# PolSAR_GWLS_filter
This is a Matlab implementation demo of PolSAR GWLS filter.  
In this demo, some operations are added to avoid computing exceptions, and the efficiency of the code is also improved.          
<br/>  
Please refer to the following paper for a more detailed description of the algorithm.  
Yexian Ren, Jie Yang, Lingli Zhao, Pingxiang Li, Zhiqu Liu, Lei Shi. A Global Weighted Least Squares Optimization Framework for Speckle Filtering of PolSAR Imagery. IEEE Trans. Geosci. Remote Sens. （https://ieeexplore.ieee.org/document/8457251）  
<br/> 
1. The **AutoParaforWLS2.m** is used to determine parameters automatically.  
<br/> 
2. The **spanWLSfilterV2.m** is used to do the global weighted least squares fitting.  
<br/> 
3. The **fPolRGBshow.m** is used to visualize PolSAR images and it is helpful to visualize high resolution PolSAR images. We have used this code in the following paper for the high resolution PolSAR image visualization.   
Ren et al. "SIRV-Based High-Resolution PolSAR Image Speckle Suppression via Dual-Domain Filtering" ,IEEE Trans. Geosci. Remote Sens.

<br/>   

**ACKNOWLEDGMENT**    
Thank Dr. Xianxiang Qin for his constructive suggestions to improve this work.       
<br/>   

Thank you for your reading！  

Yexian Ren  
Email: renyexian@foxmail.com  
QQ: 2538715345  
