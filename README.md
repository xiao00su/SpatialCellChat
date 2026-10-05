Spatial CellChat (i.e., CellChat v3) is an updated version of CellChat toolkit that

- enables the [inference of cell-cell communication at single-cell resolution from spatial transcriptomics data](https://htmlpreview.github.io/?https://github.com/jinworks/SpatialCellChat/blob/master/tutorial/SpatialCellChat_analysis_of_spatial_transcriptomics_data.html). 
- is applicable to diverse technologies of spatial transcriptomics data. We add [Frequently Asked Questions (FAQ) when analyzing spatially resolved transcriptomics datasets](https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html), particularly on how to apply Spatial CellChat to different technologies of spatial transcriptomics data, including sequencing-based and in-situ imaging-based readouts.

To ensure efficient and scalable inference of cell-cell communication at single-cell resolution from spatial transcriptomics data, Spatial CellChat optimizes the data structure within CellChat object. To enable users still can run their previously calculated CellChat v1/v2 object and smoothly upgrade to CellChat v3, we currently deposite the source codes and tutorials of Spatial CellChat at this GitHub repository. 

The instructions of installation and tutorials are available at [CellChat toolkit](https://github.com/jinworks/CellChat).

setwd('/home/dell/下载/SpatialCellChat-main')
devtools::document() # 生成文档
readLines("NAMESPACE") # 检查NAMESPACE文件内容
devtools::build() # 构建包
pak::local_install() # 本地安装测试


