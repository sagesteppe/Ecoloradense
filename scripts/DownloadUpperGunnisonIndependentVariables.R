# the domain for all analyses will be the Upper Gunnison, where 3m resolution data are 
# available. 
# From the RMBL Spatial Data Platform we will download an Digital Elevation Model
# https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release3/UG_dem_3m_v1.tif

# From this dataset we will derive several variables: 
# Aspect (WhiteBoxTools), which we will decompose into Southness and Northness  
# Slope  (WhiteBoxTools), 
# Topographic position index  (WhiteBoxTools) 
# Topographic roughness index  (WhiteBoxTools) 
# Roughness (Terra) 
# Geomorphons  (WhiteBoxTools) 
# Profile Curvature (WhiteboxTools)   
# Plan Curvature (WhiteboxTools)   
# Compound Topographic Index   

# From NAIP data we calculate the following metrics  
# Percent fractional rock/bedrock  
# Percent soil cover   
# Canopy height model 
# NDVI 

# Other products downloaded the RMBL spatial data platform include: 
# summer subcanopy solar radiation 
# https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release2/UER_srad_subcanopy_day265_3m_v1.tif
# 

