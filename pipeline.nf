input = file(params.input)
genome = file(params.genome)
//genomeIndex = file(params.genomeIndex)
//fai = file(params.fai)
/*
motifs = file(params.motifs)
motifDb = file(params.motifDb)
chromSizes = file(params.chromSizes)
genome = file(params.genome)
phastCons = file(params.phastCons)
gtf = file(params.gtf)
*/

(fastqs, fastqs2) = Channel
.from(input.readLines())
.map { line ->
  list = line.split()
  id = list[0]
  fastq = file(list[1])
  adapter = list[2]
  [ id, fastq, adapter ]
}.into(2)



process indexGenome {
	
	input:
	file genome
	
	output:
	set file("genomeDir") into genomeIndex
	set file("${genome}.fai") into fastaIndex
	
	script:
	"""
	mkdir genomeDir
	samtools faidx ${genome}
	STAR --runThreadN 18 --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles ${genome}
	"""
}

//process faidx {
//	
//	input:
//	file genome
//	
//	output:
//	set file("genomeDir") into genomeIndex
//	set file("${genome}.fai") into fastaIndex
//	
//	script:
//	"""
//	mkdir genomeDir
//	samtools faidx ${genome}
//	STAR --runThreadN 18 --runMode genomeGenerate --genomeDir genomeDir --genomeFastaFiles ${genome}
//	"""
//}

/*
process indexHairpin {

	input:
	
	output:
}
*/

//process projExons {
//	
//	input:
//	
//}

process cutadapt {

	input:
    set expId, file(fastq), adapter from fastqs

	output:
	set expId, file('trim.fastq.gz') into trimmedFastq
	
	script:
//	fastqsize=$(stat -c%s "${fastq}")
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
	set expId, file(fastqTrim) from trimmedFastq

	output:
	set expId, file('Aligned.sortedByCoord.out.bam') into bamGenome

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

(bamGenome1, bamGenome2) = bamGenome.into(2)

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

/*


process bamGenome2counts {

	input:
	file fai
    set expId, file(bam) from bamGenome2

	script:
	"""
	count.reads.elements.py --abam ${bam} -b  -o ${expId}.gene.reads.tsv
	"""
}




motifSeq = Channel
.from(motifs.readLines())
.spread(peakSeq)

process fimoScan {
	
	publishDir "${outDir}/${expId}/${motif}", mode: 'link'

	input:
	file motifDb
	set motif, expId, file(peakSeq) from motifSeq

	output:
	set expId, motif, file('fimo.txt') into fimoOut

	script:
	"""
	fimo --motif ${motif} --text --parse-genomic-coord ${motifDb} ${peakSeq} | awk '\$6>0' > fimo.txt
	"""
}

(fimoOut1, fimoOut2) = fimoOut.filter { it ->
   it[2].readLines()
   .findAll { line ->
       line =~ /^[^#]/
   }.size() > 0
}.into(2) 

motifReads = bigwigs 
.cross(fimoOut1) 
.map { bw, fimo ->
  fimo + [ bw[1], bw[2] ]
}

process motifMatrix {

	publishDir "${outDir}/${expId}/${motif}", mode: 'link'

	input:
	set expId, motif, file(fimo), file(bwFwd), file(bwRev) from motifReads 

	output:
	set expId, motif, file('read.starts.matrix'), file('fimo.txt') into (matrixStarts1, matrixStarts2)

	script:
	"""
	paste <(cat ${fimo} | tail -n+2 | cut -f2-4 | bwtool -pseudo=0 -fill=0 matrix 100:100 stdin ${bwFwd} stdout) <(cat ${fimo} | tail -n+2 | cut -f2-4 | bwtool -pseudo=0 -fill=0 matrix 100:100 stdin ${bwRev} stdout) > read.starts.matrix
	"""

}

process annMotif {

	publishDir "${outDir}/${expId}/${motif}", mode: 'link'

	input:
	file phastCons
	file gtf
	set expId, motif, file(fimo) from fimoOut2

	output:
	set expId, motif, file('ann.bed') into annoFimo

	script:
	"""
	cat ${fimo} | tail -n+2 | cut -f2-4 | bwtool summary -pseudo=0 -fill=0 stdin ${phastCons} stdout | cut -f8 > ConsScore.txt
	
	cat ${gtf} | awk '\$3==\"transcript\"' | extendBed.py --gtf --starts | sort -k 1,1 -k4,4n | closestBed -d -t first -a <(cat ${fimo} | tail -n+2 | cut -f2-4) -b stdin | awk -F \"\\t\" '{print \$NF}' > TSSdist.txt
	
	paste <(cat ${fimo} | cut -f2-6 | tail -n+2) ConsScore.txt TSSdist.txt | sed '1ichrom\thg19Start\thg19End\tStrand\tPWMscore\tConsScore\tTSSdist' | awk '\$5>=0' > ann.bed
	"""
}


process centSimple {

	publishDir "${outDir}/${expId}/${motif}", mode: 'link'

	when:
	params.model == 'simple'

	input:
	set expId, motif, file(matrix), file(fimo) from matrixStarts1

	output:
	set file('annSimple.centFit.bed'), file('PostPr.simple.txt') into centSimpleOut

	script:
	"""
	centipede.R -m ${matrix} --DampLambda 0.1 --model simple -o PostPr.simple.txt
	paste <(cat ${fimo} | cut -f2-6 | tail -n+2) PostPr.simple.txt | sed '1ichrom\thg19Start\thg19End\tStrand\tPWMscore\tPostPr' > annSimple.centFit.bed
	"""

}


matrixStartsAll = matrixStarts2
.cross(annoFimo) { it ->
  [ it[0], it[1] ]
}
.map { matrix, anno ->
  matrix + [ anno[2] ]
}

process centAll {

	publishDir "${outDir}/${expId}/${motif}", mode: 'link'

	when:
	params.model == 'all'

	input:
	set expId, motif, file(matrix), file(fimo), file(annoFimoBed) from matrixStartsAll

	output:
	set file('annAll.centFit.bed'), file('PostPr.all.txt') into centAllOut

	script:
	"""
	centipede.R -m ${matrix} --DampLambda 0.1 --model all -a ${annoFimoBed} -o PostPr.all.txt
	paste <(cat ${fimo} | cut -f2-6 | tail -n+2) PostPr.all.txt | sed '1ichrom\thg19Start\thg19End\tStrand\tPWMscore\tPostPr' > annAll.centFit.bed
	"""

}

*/

workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Execution status      : ${ workflow.success ? 'OK' : 'failed' }"
    println "Duration              : ${workflow.duration}"
}
