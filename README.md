# MinSA
A Matlab program for calculating the mean relative error (MRE) under different numbers of view-fields (fieldN) and areas of view-field (fieldS).

### _______________________________________________________________________________________________
### Research paper: 
### Effect of number and area of view-fields on the measurement accuracy of hair follicle density in goats (Capra hircus)
# Abstract
Improving the hair follicle density (HFD) has long been a breeding goal, as it can directly affect the fiber yield of fiber-bearing animals, such as cashmere goats, sheep, rabbits and alpacas. Due to the limitation of sample size and time consumption, the HFD can not be easily obtained from a large histological area. If the area is too small, the experimental sampling error cannot be ignored. The number of view-fields (fieldN), the area of view-field (fieldS) and the total statistical area (total SA, fieldN×fieldS) needed in the HFD examination have not been determined. Thus, we aimed to develop a strategy to address this unclear topic. Goat-skin was taken as materials, and paraffin sectioning and hematoxylin-eosin staining were done. Images of each skin sample were captured, stitched and clipped. All hair follicles (HFs) were manually marked as color dots, and the background noise was removed. We acquired 52 clean images (50 mm2 of stained skin section). We considered 10 types of fieldN, fieldN∈[1,10] and 14 types of fieldS, from 32.0 mm2 to 0.10 mm2. A Matlab-based algorithm/program, MinSA, was developed to calculate the mean relative error (MRE) under different fieldN and fieldS. When fieldN=1, the field area was 5.04 mm2 where MRE reached 5%, indicating that if the field area was less than 5.04 mm2 (one field) for the HFD examination, the MRE would be higher than 5%. As the fieldN increased from 2 to 10, the fieldS at 5% MRE gradually reduced, but the total SA increased from 4.45 to 5.49 mm2. The results demonstrated that using multiple view-field was not a good choice. For goats, 1-3 fields are recommended, and the total statistical area should be over 5 mm2 of stained skin section. Furthermore, this algorithm could apply to other economic fur-bearing animals.
