# nucleServ

Scripts and wrappers to run nucleR and NucleosomeDynamics on the backend of
a server

# Attributes

# Stiffness

Primary data
Score: Stiffness estimation. It's the energy required to move the nucleosome (expressed in cal/bp), it's derived from the gaussian standard deviation.

Attributes
nucleR_score: the nucleR score given to that nucleosome
nucleR.class: the nucleR class given to that nucleosome
gauss_k: the height of the peak of the gaussian curve
gauss_m: the position of the peak of the gaussian curve
gauss_sd: the standard deviation of the gaussian curve


# NucleosomeDynamics

Primary data
Position: region where the movement happens
Type: change in the nucleosome map
Score: magnitude of the change

Attributes
class: type of hotspot (see help for all possible types)
nuc: to which nucleosome the movement belongs. NA means that the hostpot couldn't be unequivocally associated to one nucleosome.
number_of_reads: number of reads involved in this movement
hreads: number of reads involved in the movement relative to the number of reads present in the area. This value ranges from 0 to 1 and the closest it is to 1, the more significant the movement.


# nucleR

Primary data
Score: Positionning score. It is calculated as the weighted sum of width and height scores.

Attributes
score_width: Witdth score. It is a measure of how sharp a peak is. A value of 0 would be an extremely wide peak and a value of 1 a very sharp one.
score_height: Height score. Tells how large a peak of a nucleosome is. The bigger this number, the higher the peak.
class: Whether the nucleosome is well-positioned (W) or fuzzy (F) or undetermined. The taken value depends on score_h and score_w. Undetermined means the exact position of the nucleosome cannot be determined due to strong fuzziness.


# Periodicity

Attributes
nucleosome_first: First nucleosome of the gene.
nucleosme_last: Last nucleosome of the gene.
score_phase: Is a measure of the phase between the first and the last nucleosome. A score of 0 means the nucleosome are completely phased and a score of 82 corresponds to totally antiphased nucleosomes.
score_autocorrelation: It is directly computed from the experimental coverage and is quantitative measure of the periodicity of nucleosomes inside the gene body.


# TSS classes

Primary data
Position: Region between two nucleosomes surrounding the TSS.

Attributes
classification: Descriptor of the Transcription Start Site. See the help for possible options.
distance: Distance in base pairs between the nucleosome +1 and the nucleosome -1.
id
nucleosome minus1: Position of the nucleosome -1.
nucleosome plus1: Position of the nucleosome +1
TTS_position: Position of the Transcription Start Site.
