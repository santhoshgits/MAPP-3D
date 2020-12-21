# MAPP-3D is a New method for aligning protein binding sites ( both pairwise and multiple site alignment ) written in Python version 2.7.
**Two types of script has been made to serve two independent purpose**
1) **pocket_matrix7.py** - This code carries out binding site alignment taking a pair of binding sites as the inputs (**i.e pairwise alignment**).
It takes atomic coordinates of the binding sites in PDB format as the input, finds optimal alignment, and superposes site1 onto site2.<br>
Usage: python2.7 pocket_matrix7.py <site1.pdb> <site2.pdb><br>
Example: python2.7 pocket_matrix7.py 1L2T_ATP_A_1301.pdb 1TF7_ATP_A_1901.pdb 

```markdown
1. align.txt - contains list of residues matched between sites
2. fixed.pdb - the translated coordinate of fixed binding site
3. frag.pdb - the translated coordinate of reference site
4. site1.pdb - same as fixed.pdb but contains only the coorinates of matched residues
5. site2.pdb - same as frag.pdb but contains only the coorinates of matched residues
```
###### All required files and scripts can be found in the directory 'PairWiseComparison'

2) **pocket_matrix7_mpi.py** - This is the extension the above script 'pocket_matrix7.py' with added MPI support. This code is designed to handle large number of comparisons (i.e. million pairs ) on a cluster or supercomputer or multi-core desktop machines.

> ### @icon-info-circle NOTE
>This module only report two things, i) residue-residue correspondence and ii) MAPP scores for all input pairs. It doesn't output the newly aligned coordinates. However, if a user is interested to see the translated coordinated of multiple alignment in a graphics viewer such as Pymol, we recommend them to use pairwise alignment code. 


###### All required files and scripts for this mode can be found in the directory 'MulitpleSiteAlignment'

---

<span style="color:blue">INSTALLATION INSTRUCTION</span>
1) To run pocket_matrix7.py, no installation is required other than python v-2.7.
2) For parallel implementation of MAPP i.e pocket_matrix7_mpi.py, user have to install MPI librarires like MPICH compatible for their system architecture.
Please refer to this link regarding installation instruction https://www.mpich.org/downloads/
After successful installation of mpich. User then have to install a python MPI package 'mpi4py' to talk with main MPI binaries.
The command for installing mpi4py is 'pip install mpi4py'.

Note: All the above softwares will require sudo previlige for installations. So please invoke terminal with sudo command.

---

---
<h3><span style="color:red">Handling ERRORS</span></h3>

1) Both version of MAPP requires binding site coordinates in .pdb format and chain identifier shoule be present for all the residues.

2) While running MPI version of MAPP. If a progam gets terminated in between due to some reason. Please don't delete the 'align_output.txt', MAPP also reads this file as a checkpoint
and run only those pairs that are not compared before.

File align_output.txt is tab delimited

---

# Steps for Running MPI version of MAPP

For the purpose of this tutorial, we have added a total of 100 ATP binding sites to a folder 'ATP'. The goal is to run MAPP(MPI) for all ATP pairs and find if any motif could be inferred.

1. **Run Pairs.py script -** This will read all sites in the given folder and create a tab separated paired entries stored in 'Pairs.txt'.  
User can provide their custom pair wise entries to avoid pairs that are not going to be similar. Eg) When two different ligand binding site are present are in same folder
then it is advisable to create a pairs for .<br>
Usage: python Pairs.py ATP <br> 
Output: PairList.txt


2. **Run pdb_res.py -** This is sort all binding sites based on the number of residue present in them.<br> Usage: python PDBSize.py ATP <br>
Output: PDBSize.txt

3. **Running MAPP**<br>
USAGE: mpirun -n 4 python <arg1> <arg2> <arg3><br>
arg1 - ATP site folder<br>
arg2 - output of Pairs.py<br>
arg3 - output of PDBSize.py

Example) mpirun -n 4 python ATP PairList.txt PDBSize.txt <br>
align_output.txt is the output file generated after running the step 3.

4. **Analysing MAPP Result**-
File 'align_output.txt' is a tab separated data containing MAPP scores and residue-residue correspondance for all combinations of site pairs as specified in file 'PairList.txt' for which the coordinate is present in the folder 'ATP'


><span style="color:red">  @icon-info-circle To find representative for the ATP binding site, any clustering algorithm can be used. Here we pick one representative based on number connections with the other ATP sites
after imposing a MAPP score cutoff of M-dist-min > 0.6 and M-dist-max > 0.4.</span>

NOTE: The nature of binding site varies from ligands to ligands. Hence it is advisable to first analyse the site network and pick the correct representative.

Run the below provided script to find the representative site<br>
python2.7 Analyse.py align_output.txt .4 7

4. **Generating Motif based on the chosen representative**<br>
Usage: Motif.py <align_output.txt> < representative-site> <No. of residue match><br>
Example: Motif.py align_output.txt 1B0U_ATP_A_301.pdb 4

Output: [WW]-x-x-[IHF]-x-[VA]-x(18)-[PSA]-[TS]-G-[SA]-G-K-[ST]-T-x(22)-[EQ]-x(78)<br>
The sequence [TS]-G-[SA]-G-K-[ST]-T is a well characterized Walker motif associated with nucleotide binding was identified correctly by MAPP. 