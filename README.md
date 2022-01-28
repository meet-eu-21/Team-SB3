# Team-SB3

Welcome to Meet-EU Team SB3 (Sorbonne University) team !

We decided to work on ***topic B : Chromosome compartments !!***

The aim of our group is to tackle chromosome compartments, first by seing what compartments can be identified within a chromosom and then between chromosoms.


Most of the pipeline is taken from Leopold Carron in files [HiCToolBox](Code/HiCtoolBox.py) and [main](Code/main.py) 


All the HiC data, compartment gold standards and gene density files are taken from this repository : http://www.lcqb.upmc.fr/meetu/dataforstudent/

# Tutorial : how to use our code to get compartments and compare it to gold standards, featuring IFB cluster !

This present tutorial aims at helping you use our code. It covers all the steps we describe in our report, from reading Hi-C matrix to evaluation of precision compared to ~~Leopold Carron~~ gold standards.

## HiC contact maps reading and compartment generation

This covers how to use [Compartments_SB3_v_Val2.py](Code/Compartments_SB3_v_Val2.py) that you can find in the main branch

### A brief summary of this script : 
The purpose of this code is to find DNA compartments by analysis of HiC contact maps and Gene Density data. The code is operational for the study of intrachromosomal data (i.e., study of a single chromosome in contact with itself) as this part was designed by Leopold Carron and his team. As a consequence, this part will not be explained here.
What will be explained here is the study of interchromosomal data (i.e., study of the contacts of two different chromosomes) as designed by our team from the above part.

### To run the code :
At the bottom end of the script, you'll find several variables to caliber before running the code :
- A string variable mode which you can caliber to "inter" or "intra". This variable allows the script to launch the corresponding study (interchromosomal for "inter" and intrachromosomal for "intra"). For each of these modes, you can caliber the following variables independently from the other mode.
-  An integer variable resolution which has to be the resolution of your HiC data.
-  A string variable HiC_fic which is the name of your HiC data file.
-   One or two string variables depending of the chosen mode. The intrachromosomal mode has the geneDen_fic variable which is the name of the Gene Density file associated to your sole studied chromosome. The interchromosomal mode has the geneDen_fic_1 and geneDen_fic_2 variables which are the names of the Gene Density files associated to your two studied chromosome. IMPORTANT : The Gene Density files have to be given in the same order than the chromosomes are arranged in the HiC data file name.
Once you have calibrated those variables, you can run the code which will execute the corresponding pipeline on the given files at the given resolution.


### Interchromosomal study :
What follows is a description of the interchromosomal study procedure. Unless noted as such (and thereby explained), all steps of the procedure are the same as in the intrachromosomal study : 

- After loading the HiC data in an integer matrix format, this matrix is directly reformatted into a coordinate matrix, without being concatenated to its own transpose.
The reason for this is that contrarily to the intrachromosomal HiC data, where all coordinated contacts are symmetrical (hence why we have to explicitly tell the matrix that for every couple of coordinates [x ; y], the complementary couple [y ; x] has the same value), the interchromosomal HiC data already contain complementary couples and of possible different value at that (couples [xA ; yB] and [xB ; yA] are mathematically the same but biologically different as each x and y come from a different chromosome).

- After being transformed into a binned map, the data are filtered by a modified version of the filtering function from the HiCtoolbox script. Like its parent, the modified function returns the binned map rid of its empty lines and columns (meaning there is no contact between the chromosomes for a particular coordinate) but since the studied chromosomes are different, we had to make sure the filtering isn't symmetrical (i.e., chromosome A may have contacts with chromosome B with its coordinate xA while for the same coordinate xB on chromosome B, there might be no contact with chromosome A). For the same reason, the modified function returns as well the saved bins lists of both chromosomes.
The modified version of the filtering function is findable at the top of the script, above the intrachromosomal pipeline.

- Next is the applying of the Sequential Component Normalization (SCN for short). The function called is once again a modified version of the one from the HiCtoolbox script. The modification here consists in returning the normalized matrix as it is, without multiplying it with its transpose.
The reason for this is to keep the matrix non squared. Indeed, as the studied chromosomes may not be of the same length, the matrix's dimensions may not be equal from the start, so we have to keep it this way..
The modified version of the SCN function is findable at the top of the script, above the modified filtering function.

At this point of the study, the pipeline generates its first visualized result, the initial heatmap.

- Then, we compute the correlation for Observed over Expected contacts (O/E matrix for short). For the same reason as for the SCN, the O/E computation is modified to take into account the non-squareness of the matrix by use of both dimensions.
The visualization of the O/E heatmap is then generated.

Next is the computation of the correlation matrices, plural. Indeed, because we are working with two different chromosomes, we need to correlate the data from each chromosome's point of view. As a result, we find ourself with two visualized correlation heatmaps.

- As a consequence of having two correlation matrices, all of the remaining steps are doubled : the eigenvalues and eigenvectors calculation, the Gene Density data analysis (coupled with the use of both list of saved bins from the filtering) and the finding of the compartments.
However, the data processing for all of this stay the same as in the intrachromosomal study.

Our adaptation of the intrachromosomal study to the interchromosomal case ends here.
Technically speaking, the original code has a last step in which it generates a .pdb file. This step has not been adapted yet in our code, that is why it is overlooked during the code's execution by being between quotation marks. But the procedure is the same as it is for intrachromosomal study. 

## Compartments evaluation to gold standards :

This part covers how to use [compare results][Code/compare_results.py] to convert the file of compartments generated by the previous pipeline to get similarity compared to gold standard results given by Leopold in the folder [Gold standards][Gold_standards/]. 

### The general structure of the script :

It is composed of two main methods, one for intrachromosomal and the other for interchromosomal comparision with gold standards.

### How to launch the code :

Just load the file with your favourite IDE and then type :


## Collaborative work with team SB1 (Sorbonne University) : 

## Launching our code on the IFB cluster : 

## License

This code is distributed under the MIT license. See [LICENSE][LICENSE.txt] for more information.
