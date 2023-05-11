# scATAC_code
scATAC_code workflow


1. 《code of 2021cell_stem_cell》未能绕过snapATAC，中间文件出图  
2. 《code of 2022nature_immu》无代码，使用该文章的原始数据  
3. 《code of 2022scitfic_report》虽然本流程可用，但docker太慢，且过程中容易出小问题，使用ArchR可完全重复文章图片  
4. 《code of 2023cancer_res》本文提供了完整代码，但是数据为bam文件，接下来使用nature immu的数据复现。现只完成了ArchR对数据的筛选处理，并未走Signac步骤，后续可继续进行（ArchR卡在call peak步骤）