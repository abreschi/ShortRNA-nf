input = file(params.input)
genome = file(params.genome)
annotation = file(params.annotation)


// Read input file
fastqs = Channel
.from(input.readLines())
.map { line ->
  list = line.split()
  id = list[0]
  fastq = file(list[1])
  adapter = list[2]
  [ id, fastq, adapter ]
}


/*
	Common processes
*/

process indexGenome {
	
	input:
	file genome
	
	output:
	file "genomeDir" into genomeIndex
	
	script:
	"""
	mkdir genomeDir
	STAR --runThreadN 18 --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles ${genome}
	"""
}

process faidx {
	
	input:
	file genome
	
	output:
	file "${genome}.fai" into (fastaIndex, fastaIndex0, fastaIndex1)
	
	script:
	"""
	samtools faidx ${genome}
	"""
}


process indexHairpin {

	input:
	file fai from fastaIndex0
	file genome
	file annotation
	
	output:
	file "genomeDir" into hairpinIndex
	file "hairpins.gtf" into (hairpinsGtf1, hairpinsGtf2)
	
	script:
	"""
	mkdir genomeDir

	cat ${annotation} | awk '\$3==\"gene\"' | grep "gene_type \\\"miRNA\\\"" > hairpins.gtf

	cat hairpins.gtf | awk 'BEGIN{OFS=FS=\"\t\"}{n=split(\$9,a,\"; \"); \
	for(i=1;i<=n;i++){split(a[i],b,\" \"); gsub(/\"/, \"\", b[2]); \
	if(b[1]==\"gene_id\"){gn=b[2]}} print \$1, \$4, \$5, gn, 0, \$7}' \
	| fastaFromBed -fi ${genome} -bed stdin -name -s -fo hairpin.fa

	samtools faidx hairpin.fa

	Nbases=\$(cat hairpin.fa.fai | awk '{c+=\$2}END{a=((log(c)/log(2))/2 -1); print (a>14?14:a)}')

	STAR --runThreadN 16 --runMode genomeGenerate --genomeSAindexNbases \$Nbases --genomeDir genomeDir --genomeFastaFiles hairpin.fa
	"""
}


process projExons {
	
	input:
	file annotation
	
	output:
//	file annotation.baseName.replace(".gtf", ".exonproj.gff")
	file "exonproj.gff" into exons

	script:
	"""
	awk '\$3==\"exon\"{gn=\$10; tr=\$12; \$10=tr; \$12=gn; \$9=\"transcript_id\"; \$11=\"gene_id\"; print}' ${annotation} | cutgff.awk | gff2gff.awk > formakeSP.gff
	makeSP formakeSP.gff -f gff -v > formakeSP.gff_segproj.gff
	awk '{split(\$10,a,\":\"); print a[1], \"makeSP\", \"projex\", a[2], a[3], \".\", (a[4]==\"p\" ? \"+\" : \"-\"), \".\", \"gene_id\", \"\\""\$12\"\\";\";}' formakeSP.gff_segproj.gff | sort | uniq | gff2gff.awk > exonproj.gff
	"""
	
}


/* 
	Individual file
*/

process cutadapt {

	time { fastq.size() > 20_000_000_000 ? 12.h : 6.h }

	input:
    set expId, file(fastq), adapter from fastqs

	output:
	set expId, file('trim.fastq.gz') into (trimmedFastq1, trimmedFastq2)
	
	script:

	"""
	if [ ${adapter} != 'polyA' ]; then
		cutadapt -a $adapter --trim-n -o trim.fastq.gz -m 16 --too-short-output too-short.fastq.gz $fastq 
	else
		cutadapt -a \"A{10}\" -n 10 -e 0.1 --trim-n -o trim.fastq.gz --too-short-output too-short.fastq.gz $fastq
	fi
	"""
}


process map2genome {

	input:
	set file(genomeDir) from genomeIndex.first()
	set expId, file(fastqTrim) from trimmedFastq1

	output:
	set expId, file('Aligned.sortedByCoord.out.bam') into (bamGenome1, bamGenome2)

	script:
	"""
	params='--readFilesCommand zcat \
		--outReadsUnmapped Fastx \
		--outFilterMultimapNmax 10 \
		--outFilterMultimapScoreRange 0 \
		--outSAMunmapped None \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterScoreMinOverLread 0 \
		--outFilterMatchNminOverLread 0 \
		--outFilterMatchNmin 16 \
		--outFilterMismatchNmax 1 \
		--alignSJDBoverhangMin 1000 \
		--alignIntronMax 1 \
		--runThreadN 16 \
		--limitBAMsortRAM 30000000000
	'
	STAR --genomeDir ${genomeDir} --readFilesIn ${fastqTrim} \$params
	samtools index Aligned.sortedByCoord.out.bam
	"""
}


process bamGenome2bw {

	input:
	set file(fai) from fastaIndex.first()
    set expId, file(bam) from bamGenome1

	script:
	"""
	genomeCoverageBed -split -bg -strand '+' -ibam ${bam} > ${expId}.plusRaw.bg
	bedGraphToBigWig ${expId}.plusRaw.bg ${fai} ${expId}.plusRaw.bw
	rm ${expId}.plusRaw.bg

	genomeCoverageBed -split -bg -strand '-' -ibam ${bam} > ${expId}.minusRaw.bg
	bedGraphToBigWig ${expId}.minusRaw.bg ${fai} ${expId}.minusRaw.bw
	rm ${expId}.minusRaw.bg
	"""
}



process bamGenome2counts {

	input:
	file prExons from exons.first()
    set expId, file(bam) from bamGenome2

	script:
	"""
	count.reads.elements.py --abam ${bam} -b ${prExons} -o ${expId}.gene.reads.tsv
	"""
}

process map2hairpins {

	input:
	file genomeDir from hairpinIndex.first()
	set expId, file(fastqTrim) from trimmedFastq2

	output:
	set expId, file('Aligned.sortedByCoord.out.bam') into bamHairpin

	script:
	"""
	params='--readFilesCommand zcat \
		--outReadsUnmapped Fastx \
		--outFilterMultimapNmax 10 \
		--outFilterMultimapScoreRange 0 \
		--outSAMunmapped None \
		--outSAMtype BAM SortedByCoordinate \
		--outFilterScoreMinOverLread 0 \
		--outFilterMatchNminOverLread 0 \
		--outFilterMatchNmin 16 \
		--outFilterMismatchNmax 1 \
		--alignSJDBoverhangMin 1000 \
		--alignIntronMax 1 \
		--runThreadN 8 \
		--limitBAMsortRAM 30000000000
	'
	STAR --genomeDir ${genomeDir} --readFilesIn ${fastqTrim} \$params
	samtools index Aligned.sortedByCoord.out.bam
	"""
}


process hairpins2genome {
	
	input:
	set expId, file(bam) from bamHairpin
	file gtf from hairpinsGtf1.first()
	file fai from fastaIndex1.first()

	output:
	set expId, file("${expId}.gCoord.bam") into bamHairpinGenome

	script:
	"""
	hairpin2genome.sh ${gtf} ${bam} ${fai} ${expId}.gCoord
	"""
}

process hairpins2counts {

	input:
	set expId, file(bam) from bamHairpinGenome
	file gtf from hairpinsGtf2.first()

	script:
	"""
	count.reads.elements.py --abam ${bam} -b ${gtf} -o ${expId}.reads.tsv
	"""
}

workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Execution status      : ${ workflow.success ? 'OK' : 'failed' }"
    println "Duration              : ${workflow.duration}"
}
