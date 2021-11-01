# SiteMotif (MAPP-3D) is a graph based method for aligning protein binding sites ( both pairwise and multiple site alignment ) in sequence-order independent fashion written in Python version 2.7.
**Two types of scripts have been made to serve two independent purposes**
1) **pocket_matrix7.py** - This code perform binding site alignment taking a pair of binding sites as the inputs (**i.e pairwise alignment**).
It takes input an atomic coordinates of the binding sites in PDB format, finds optimal alignment, and superposes site1 onto site2.<br>
Usage: python2.7 pocket_matrix7.py <site1.pdb> <site2.pdb><br>
Example: python2.7 pocket_matrix7.py 1L2T_ATP_A_1301.pdb 1TF7_ATP_A_1901.pdb 

```markdown
1. align.txt - contains list of residues matched between sites
2. fixed.pdb - translated coordinate of the fixed binding site
3. frag.pdb - translated coordinate of the reference site
4. site1.pdb - same as fixed.pdb but contains only the coorinates of matched residues
5. site2.pdb - same as frag.pdb but contains only the coorinates of matched residues
```
###### All required files and scripts can be found in the directory 'PairWiseComparison'

2) **pocket_matrix7_mpi.py** - This is the extension the above script 'pocket_matrix7.py' with an added MPI support. This code is designed to handle large number of comparisons (i.e. million pairs ) on a cluster or supercomputer or multi-core desktop machines.

> ### NOTE
>This module only report two things, i) residue-residue correspondence and ii) M-scores for all input pairs. It doesn't output the newly aligned coordinates. However, if a user is interested to see the translated coordinated of multiple alignment in a graphics viewer such as Pymol, we recommend them to use pairwise alignment code. 


###### All required files and scripts for this mode can be found in the directory 'MulitpleSiteAlignment'

---

<span style="color:blue">INSTALLATION INSTRUCTION</span>
1) To run pocket_matrix7.py, no installation is required other than python version 2.7.
2) For parallel implementation of SiteMotif i.e pocket_matrix7_mpi.py. User have to install MPI librarires such as MPICH compatible for their system architecture.
Please refer to this link regarding installation instruction https://www.mpich.org/downloads/.
After successful installation of mpich. User then have to install a python MPI package 'mpi4py' to communicate with main MPI binaries.
The command for installing mpi4py is 'pip install mpi4py'.

Note: All the above softwares will require sudo previlige for installations. So please invoke the terminal with sudo command.

---

---
<h3><span style="color:red">Handling ERRORS</span></h3>

1) Both version of SiteMotif requires binding site coordinates in .pdb format and chain identifier have to be present for all the residues.

2) While running MPI version of our code. If a progam gets terminated in between due to some unknown reason. Please don't delete 'align_output.txt', pocket_matrix7_mpi.py read  file 'align_output.txt' as a checkpoint and runs only those pairs that are not compared before.

File align_output.txt is tab delimited

---

# Steps for Running MPI version of MAPP

For the purpose of this tutorial, we have added a total of 100 ATP binding sites to a folder 'ATP'. The goal is to run SiteMotif(MPI) for all ATP pairs and find if any motif could be inferred.

1. **Run Pairs.py script -** This will read all sites in the given folder and create a tab separated paired entries stored in 'Pairs.txt'.  
User can provide their custom pair wise entries to avoid pairs that you dont want to compare.
Usage: python2.7 Pairs.py ATP <br> 
Output: PairList.txt
PairList.txt is tab separated


2. **Run pdb_res.py -** This is sort all binding sites based on the number of residue present in them.<br> Usage: python PDBSize.py ATP <br>
Output: PDBSize.txt

3. **Running SiteMotif**<br>
USAGE: mpirun -n 4 python2.7 pocket_matrix_mpi7.py arg1 arg2 arg3<br>
arg1 - ATP site folder<br>
arg2 - output of Pairs.py<br>
arg3 - output of PDBSize.py

Example) mpirun -n 4 python2.7 pocket_matrix_mpi7.py ATP PairList.txt PDBSize.txt <br>
align_output.txt is the output file generated after running the step 3.

4. **Analysing SiteMotif Result**-
File 'align_output.txt' is a tab separated data containing MAPP scores and residue-residue correspondances for all combinations of site pairs as specified in file 'PairList.txt' for which the coordinate is present in the folder 'ATP'


><span style="color:red"> To find representative for the ATP binding site, any clustering algorithm can be used. Here we pick one representative based on number connections with the other ATP sites
after imposing a cutoff of M-dist-min > 0.6 and M-dist-max > 0.4.</span>

NOTE: The nature of binding site varies from one ligand to the next, as well as from one binding sites to another. As a result, it is recommended to do an initial analysis of the site's network and select an acceptable representative.

Run the below provided script to find the representative site<br>
python2.7 Analyse.py align_output.txt .4 7

4. **Generating Motif based on the chosen representative (For undermining purposes only)**<br>
Usage: Motif.py <align_output.txt> < representative-site> <No. of residue match><br>
Example: python2.7 Motif.py align_output.txt 1B0U_ATP_A_301.pdb 4

Output: [WW]-x-x-[IHF]-x-[VA]-x(18)-[PSA]-[TS]-G-[SA]-G-K-[ST]-T-x(22)-[EQ]-x(78)<br>
The sequence [TS]-G-[SA]-G-K-[ST]-T is a well characterized Walker motif associated with nucleotide binding was identified correctly by SiteMotif. 
