# SoFoCompBio22 - work in progress
If you miraculously came here from my physical copy of my research report, hi! Hope you enjoyed it and hopefully more progress will have been made with my project.

If you stumbled upon this repo you may be confused. SoFoCompBio is abbreviation for "SOmmarFOrskarskola i berakningsbiologi och bioinformatik" which is a computational biology research internship at karolinska institutet, and I had the pleasure of attending its first iteration in the summer of 2022 among 5 students picked nationally in Sweden, at the Karlsson Hedestam/Murrell Lab. I was additionally even more fortunate to be accepted to the "sommarforskaskola med biomedicinskt inriktning, " karolinska institutets's long-running wet lab summer research program. This allowed me the exclusive opportunity to conduct a joint wet-lab & computational research project.

The project was about progressing the camelid germline VHH repertoire through NGS of PBMC mRNA samples from a huarizo and a computational genome mining algorithm.

The genome mining algorithm has been uploaded as a preliminary package at https://github.com/Qile0317/KmerGMA.jl

# Project abstract - (Project is unfinished)

The alpaca adaptive immune system partially produce heavy
chain only antibodies, characterized by variable and constant
regions referred to as VHH and CHH. Procurement of a com-
prehensive and diverse alpaca VHH gene repertoire is essential
for the understanding of B cell biology and has numerous ad-
vantages and benefits in therapeutics and research via usage
of alpaca nanobodies. However, the full repertoire is far from
complete. Here, we contribute to the repertoire via the creation
and execution of a modified 5’RACE protocol for next genera-
tion sequencing of both VHH and conventional VH mRNA tran-
scripts from a huarizo (Vicugna pacos × lama glama). The re-
sulting sequenced repertoire revealed over 600 thousand high
quality VDJ & VHH transcripts, including 300 thousand IgM
transcripts that can be processed in subsequent studies with
germline inference tools and experimental verification. Rudi-
mentary phylogenetic and V-gene assignment analyses pointed
strongly at the existence of novel germline alleles in our se-
quenced repertoire compared to the IMGT databse, and sub-
sequent analyses suggested a strong correlation of isotype fre-
quency to presence of nanobody hallmark mutations. Addition-
ally, we propose a novel, swift genome mining algorithm for V
gene discovery. Our unoptimized Julia implementation of the
algorithm was applied over camelid V gene loci and the Vic-
Pac3.2 full coverage alpaca genome and successfully found both
exact and approximal matches in linear time. Conclusively, our
study utilized 2 approaches to successfully progress the current
alpaca V gene repertoire to completion with high improvement
potential

# Overview of the Repo

- Data： some of the sequence data used in the paper that werent too large in filesize. 
- Figures: a collection of cool figures generated in the analysis that didnt make it to the paper due to the word limit.
- AnalysisScripts: collection of scripts (except Vsearch commands) described in my methods section that used to process the rep-seq data.
