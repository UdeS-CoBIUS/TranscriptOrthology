
# :dna: Inferring clusters of orthologous and paralogous transcripts

***Algorithm to infer clusters of isoorthologous transcripts using gene-level homology relationships and a Reciprocal Best Hits approach***


<p align="center">
<img src='./theme.png' alt='theme' width=auto height=300>
</p>

:busts_in_silhouette: __Authors__
*Wend Yam Donald Davy Ouedraogo & Aida Ouangraoua, CoBIUS LAB, Department of Computer Science, Faculty of Science USherbrooke,  Sherbrooke, Canada*

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

*   __python3 (at leat python 3.6)__
*   __NetworkX__
*   __Pandas__
*   __Numpy__
*   __ETE toolkit__


![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)

<!-- Getting started -->
<h2 id="getting-started"> :rocket: Getting Started</h2>


<!-- Main Command -->
<h3 id="main"> :computer: Main Command</h3>

> Command

<pre><code>$ python3 transcriptOrthology.py [-talg transcripts alignment] [-galg genes alignment] [-gtot gene to transcripts mappings] [-nhxt NHX gene tree] [-lowb lower bound] [-outf output folder]</code></pre>

>> Details

<table>
<tr>
    <th>argument</th>
    <th>definition</th>
    <th>format</th>
  </tr>
  <tr>
    <td>-talg | --transcriptsalignment</td>
    <td>MSA of transcripts</td>
    <td>FASTA format</td>
  </tr>
  <tr>
    <td>-galg | --genesalignment</td>
    <td>MSA of genes</td>
    <td>FASTA format</td>
  </tr>
  <tr>
    <td>-gtot | --genetotranscripts</td>
    <td>mappings g(t)</td>
    <td>FASTA format >id_transcript:id_gene</td>
  </tr>
  <tr>
    <td>-nhxt | --nhxtgenetree</td>
    <td>gene tree</td>
    <td>NHX format</td>
  </tr>
  <tr>
    <td>-lowb | --lowerbound</td>
    <td>a lower bound to select RBHs transcripts. By default equals to 0.5</td>
    <td>float number between 0 and 1</td>
  </tr>
  <tr>
    <td>-outf | --outputfolder</td>
    <td>folder to save results. The lauching program folder is set by default.</td>
    <td>String</td>
  </tr>
</table>

> Usage

<pre><code>python3 ./scripts/transcriptOrthology.py -talg ./execution/inputs/transcripts_alignments/ENSGT00390000000715.alg -galg ./execution/inputs/genes_alignments/ENSGT00390000000715.alg -gtot ./execution/inputs/transcripts_alignments/ENSGT00390000000715.fasta 
-nhxt ./execution/inputs/NHX_trees/ENSGT00390000000715.nwk -lowb 0.7 -outf ./execution/outputs/ </code></pre>

> 


![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/rainbow.png)

<h2 id="project-files-description"> :file_folder: Project Files Description</h2>


<h3 id="project-files-description-inputs"> :keyboard: Inputs description </h3>

__Inputs description__

> <em>tsm-computing()</em> :arrow_right: returns the similarity matrix tsm+ scores for all pairs of homologous transcripts.

> <em>t-clustering()</em> :arrow_right: returns the orthology graph of transcripts.

> <em>inferring-homologies()</em> :arrow_right: returns for each pair of homologous transcripts, their homology relationship type (recent-paralogs, ortho-paralogs or ortho-orthologs).

--- 

<h3 id="project-files-description-outputs"> :minidisc: Outputs description </h3>

__Outputs description__

> similarity matrix score that present the tsm+ score between each pair of homologous transcripts.

> orthology graph at the start of the algorithm showing only the pair relationships between recent-paralogs. (:warning:only retrieved if the number of transcripts is not greater than 20)

> orthology graph at the end of the algorithm showing all the relationships between pairs of isoorthologues.(:warning:only retrieved if the number of transcripts is not greater than 20)

> csv file resuming the information of the isoorthology-clustering.

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




