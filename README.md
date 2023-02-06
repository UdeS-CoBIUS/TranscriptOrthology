
# :dna: Inferring clusters of orthologous and paralogous transcripts

***Algorithm to infer clusters of isoorthologous transcripts using gene-level homology relationships and a Reciprocal Best Hits approach***


<p align="center">
<img src='./theme.png' alt='theme' width=auto height=300>
</p>

:busts_in_silhouette: __Authors__
`Wend Yam Donald Davy Ouedraogo & Aida Ouangraoua`, CoBIUS LAB, Department of Computer Science, Faculty of Science, Université de Sherbrooke, Sherbrooke, Canada*

> :bulb: If you are using our algorithm in your research, please cite our recent paper: __Upcoming__ 

> :e-mail: Contact: wend.yam.donald.davy.ouedraogo@usherbrooke.ca


<!-- TABLE OF CONTENTS -->
<h2 id="table-of-contents"> :book: Table of Contents</h2>

<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#about-the-project"> ➤ About The Project</a>
    <ol>
    <a href="#overview"> ➤ Overview</a></ol>
    <ol>
    <a href="#requirements"> ➤ Requirements</a>
    </ol>
    </li>
    <li><a href="#getting-started"> ➤ Getting Started</a>
    <ol><a href="#main"> ➤ Main command/Execution</a></ol></li>
    <li><a href="#project-files-description"> ➤ Project files descriptions </a>
    <ol>
    <a href="#project-files-description-inputs"> ➤ Inputs description</a>
    </ol>
    <ol><a href="#project-files-description-outputs"> ➤ Outputs description</a></ol>
    <ol><a href="#project-files-description-data"> ➤ Dataset</a></ol>
    </ol></li>
  </ol>
</details>

![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)

<!-- ABOUT THE PROJECT -->
<h2 id="about-the-project"> :pencil: About The Project</h2>


<!-- OVERVIEW -->
<h3 id="overview"> :cloud: Overview</h3>

*A graph-based method to infer isoorthology & recent paralogy relations in a set of homologous transcripts:dna:*

---


<!-- Requirements -->
<h3 id="requirements"> :hammer_and_pick: Requirements</h3>

*   __`python3 (at least python 3.6)`__
*   __`NetworkX`__
*   __`Pandas`__
*   __`Numpy`__
*   __`ETE toolkit`__


![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)

<!-- Getting started -->
<h2 id="getting-started"> :rocket: Getting Started</h2>


<!-- Main Command -->
<h3 id="main"> :computer: Main Command</h3>

> Command

<pre><code>usage: transcriptOrthology.py [-h] [-talg TRALIGNMENT] [-galg GENESALIGNMENT]
                              [-gtot GENETOTRANSCRIPTS] [-nhxt NHXGENETREE]
                              [-lowb LOWERBOUND] [-tsm TSMUNCORRECTED]
                              [-outf OUTPUTFOLDER]

parsor program parameter

optional arguments:
  -h, --help            show this help message and exit
  -talg TRALIGNMENT, --tralignment TRALIGNMENT
  -galg GENESALIGNMENT, --genesalignment GENESALIGNMENT
  -gtot GENETOTRANSCRIPTS, --genetotranscripts GENETOTRANSCRIPTS
  -nhxt NHXGENETREE, --nhxgenetree NHXGENETREE
  -lowb LOWERBOUND, --lowerbound LOWERBOUND
  -tsm TSMUNCORRECTED, --tsmuncorrected TSMUNCORRECTED
  -outf OUTPUTFOLDER, --outputfolder OUTPUTFOLDER</code></pre>

>> Details

<table>
<tr>
    <th>parameter</th>
    <th>definition</th>
    <th>value format</th>
  </tr>
  <tr>
    <td>-talg</br>--tralignment</td>
    <td>MSA of transcripts</td>
    <td>FASTA format </br>>{id_transcript}\n{sequence}</td>
  </tr>
  <tr>
    <td>-galg</br>--genesalignment</td>
    <td>MSA of genes</td>
    <td>FASTA format</br>>{id_gene}\n{sequence}</td>
  </tr>
  <tr>
    <td>-gtot </br>--genetotranscripts</td>
    <td>mappings g(t)</td>
    <td>FASTA format</br>>{id_transcript}:{id_gene}\n</td>
  </tr>
  <tr>
    <td>-nhxt</br> --nhxtgenetree</td>
    <td>gene tree</td>
    <td>NHX format</td>
  </tr>
  <tr>
    <td>-lowb </br>--lowerbound</td>
    <td>a lower bound to select RBHs transcripts. By default, equals to 0.5</td>
    <td>float between 0 and 1</td>
  </tr>
  <tr>
    <td>-outf </br> --outputfolder</td>
    <td>folder to save results. The current program folder is set by default.</td>
    <td>String</td>
  </tr>
</table>

> Usage example

<pre><code>python3 ./scripts/transcriptOrthology.py -talg ./execution/inputs/transcripts_alignments/ENSGT00390000000104.alg -galg ./execution/inputs/genes_alignments/ENSGT00390000000104.alg -gtot ./execution/inputs/mapping_gene_to_transcripts/ENSGT00390000000104.fasta -nhxt ./execution/inputs/NHX_trees/ENSGT00390000000104.nwk -lowb 0.7 -outf ./execution/outputs/</code></pre>

> Output expected

<pre><code>++++++++++++++++Starting ....
+++++++ All data were retrieved & the representation of subtranscribed sequences of genes into blocks are available.
+++++ Computing matrix ...       in progress
+++++ Computing matrix ...       status: Finished without errors in 1.430955410003662 seconds
+++++ Searching for recent-paralogs ...         status: processing
+++++ Searching for recent-paralogs ...         status: finished in 0.10488486289978027 seconds
+++++ Searching for RBHs ...    status: processing
+++++ Searching for RBHs ...    status: finished in 0.10737109184265137 seconds
+++++ Construction of the orthology graph (Adding nodes ...) ...        status: processing
+++++ Construction of the orthology graph (Adding nodes ...) ...        status: finished in 0.536250114440918 seconds
+++++ Searching for connected components ...    status: processing
+++++ Searching for connected components ...    status: finished in 0.536250114440918 seconds
++++++++++++++++Finished 

 Successful. Data results can be found in ./execution/outputs/
</code></pre>


![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)

<h2 id="project-files-description"> :file_folder: Project Files Description</h2>


<h3 id="project-files-description-inputs"> :keyboard: Inputs description </h3>

__Inputs files__
<ul>
    <li>:one:
    <em>tsmcomputing()</em> :arrow_right: returns the similarity matrix (tsm+ | tsm) scores for all pairs of homologous transcripts.
    <pre><code>usage: tsmComputing.py [-h] [-talg TRALIGNMENT] [-galg GENESALIGNMENT]
                       [-gtot GENETOTRANSCRIPTS] [-tsm TSMUNCORRECTED]
                       [-outf OUTPUTFOLDER]

parsor program parameter

optional arguments:
  -h, --help            show this help message and exit
  -talg TRALIGNMENT, --tralignment TRALIGNMENT
  -galg GENESALIGNMENT, --genesalignment GENESALIGNMENT
  -gtot GENETOTRANSCRIPTS, --genetotranscripts GENETOTRANSCRIPTS
  -tsm TSMUNCORRECTED, --tsmuncorrected TSMUNCORRECTED
  -outf OUTPUTFOLDER, --outputfolder OUTPUTFOLDER</code></pre>
  </li>
    <li>:two:
    <em>Tclustering()</em> :arrow_right: returns the orthology graph of transcripts.
<pre><code>usage: Tclustering.py [-h] [-m MATRIX] [-gtot GENETOTRANSCRIPTS]
                      [-nhxt NHXGENETREE] [-lowb LOWERBOUND]
                      [-outf OUTPUTFOLDER]

parsor program parameter

optional arguments:
  -h, --help            show this help message and exit
  -m MATRIX, --matrix MATRIX
  -gtot GENETOTRANSCRIPTS, --genetotranscripts GENETOTRANSCRIPTS
  -nhxt NHXGENETREE, --nhxgenetree NHXGENETREE
  -lowb LOWERBOUND, --lowerbound LOWERBOUND
  -outf OUTPUTFOLDER, --outputfolder OUTPUTFOLDER </code></pre>
  </li>
  <li>:three:
    <em>transcriptOthology()</em> :arrow_right: returns for each pair of homologous transcripts, their homology relationship type (recent-paralogs, ortho-paralogs or ortho-orthologs).
  </li>

</ul>

--- 

<h3 id="project-files-description-outputs"> :minidisc: Outputs description </h3>

__Outputs files__

<ul>
    <li>:one: <b>matrix.csv</b> : similarity matrix score that shows the tsm+ score between each pair of homologous transcripts.</li>
    <li>:two: <b>blocks_transcripts.csv|blocks_genes.csv</b> : csv file describing the representation of blocks for each transcript(resp. gene).</li>
    <li>:three: <b>start_orthology_graph.pdf|end_orthology_graph.pdf</b> : orthology graph at the start of the algorithm(resp. at the end of the algorithm) showing only the pair relationships between recent-paralogs(resp. all the orthologous clusters). (:warning:only retrieved if the number of transcripts is not greater than 20)</li>
    <li>:four: <b>groupsOfOrthoogs.csv|relationsOrthology.csv</b> : csv files resuming the information of the isoorthology-clustering.</li>
</ul>

<h3 id="project-files-description-data"> :heavy_check_mark: Dataset </h3>

***The folder data contains dataset used for the studies and also the results obtained.***

<br>
---
<br>
<br>
<br>
<br>
<br>
<br>
Copyright © 2023 CoBIUS LAB




