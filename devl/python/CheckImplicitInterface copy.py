#!/usr/bin/env python3
# -*- coding:utf-8 -*-
###
# File: \CheckFuncTypes.py
# Created Date: Wednesday, May 25th 2022, 11:51:33 am
# Author: Christian Perwass (CR/AEC5)
# <LICENSE id="Apache-2.0">
#
#   Image-Render Base Functions module
#   Copyright 2022 Robert Bosch GmbH and its subsidiaries
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
# </LICENSE>
###

class CStringAutoEnum:
    def __repr__(self) -> str:
        return f"{self.name}"

    def __str__(self) -> str:
        return f"{self.name}"

    # def __init__(self, f32_rawValue) -> str:
    #     return f"{self.name}"

    def __call__(self, f32_rawValue) -> str:
        return f"{self.name}"



if __name__ == "__main__":
    dicArgs = {
        "value": "blender/animate/object/rotate/const:1",
        "type": "label",
        "hint": "entry point identification",
        "sUnit": "entry point identification",
    }

    se =  CStringAutoEnum( "34",8 )
    print("das War's")