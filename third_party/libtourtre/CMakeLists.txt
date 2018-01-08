# PKUVIS CONFIDENTIAL
# ___________________
# 
# Copyright (c) 2009-2012, PKU Visualization and Visual Analytics Group 
# Produced at Peking University, Beijing, China.
# All rights reserved.
#                                                                              
# NOTICE: THIS IS UNPUBLISHED PROPRIETARY SOURCE CODE OF  VISUALIZATION 
# AND VISUAL ANALYTICS GROUP (PKUVIS), PEKING UNIVERSITY. DISSEMINATION
# OF  THIS  INFORMATION  OR  REPRODUCTION OF THIS  MATERIAL IS STRICTLY 
# FORBIDDEN UNLESS PRIOR WRITTEN PERMISSION IS OBTAINED FROM PKUVIS.

set (tourtre_sources
        src/ctArc.c  
        src/ctBranch.c  
        src/ctComponent.c  
        src/ctNode.c
        src/ctNodeMap.c  
        src/ctQueue.c  
        src/tourtre.c
    )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)

ADD_LIBRARY(tourtre ${tourtre_sources})

