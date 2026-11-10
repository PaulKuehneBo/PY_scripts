import os
import glob
import datetime
import random
from PIL import Image
from pathlib import Path

# <image left = "-83"  top = "-186"
#        scale = "0.639674"
#        filename = "IMG3.jpeg"
#   safeContainerLocation = "persistent" useABK = "1"> <relationships/> < /image >
GImageTemplate_str = '<image left top  filename  scale  safeContainerLocation = "persistent" useABK = "1"  > <relationships/> < /image >'


GHeightHorizontal = 550
GWidthHorizontal = 1177

# '<area areaAbilities="66299"
#       top="591.7207249381395" left="3243.228783521701"
#       height="83.27927506186052" width="706.7712164782988"
#       alpha="1" rotation="0" shadowIntensity="128" sizeborder="-1" shadowAngle="135" shadowEnabled="0" shadowBlurNew="5" areatype="textarea" zposition="1000" shadowWidthInMM="4" shadowDistance="50" colorborder="#000000"  borderenabled="0">

GTextAreaTemplate_str = '<area areaAbilities="66299" left top height width alpha="1" rotation="0" shadowIntensity="128" sizeborder="-1" shadowAngle="135" shadowEnabled="0" shadowBlurNew="5" areatype="textarea" zposition="1000" shadowWidthInMM="4" shadowDistance="50" colorborder="#000000" borderenabled="0">'
GTextAreaTemplate_str += '<text isSpineText="0"><![CDATA[<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0//EN" "http://www.w3.org/TR/REC-html40/strict.dtd"><html><head><meta name="qrichtext" content="1" /><style type="text/css">p, li { white-space: pre-wrap; }</style></head><body style=" font-family:\'Bahnschrift SemiBold\'; font-size:16pt; font-weight:400; font-style:normal;">'
GTextAreaTemplate_str += '<table style="-qt-table-type: root; margin-top:5px; margin-bottom:5px; margin-left:5px; margin-right:5px;"><tr><td style="border: none;"><p style=" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;"><span style=" color:#ffffff;">'
GTextAreaTemplate_str += (
    "$Datum$"
    + '</span></p></td></tr></table></body></html>]]><outline enabled="0" color="#000000" width="1"/>'
)
GTextAreaTemplate_str += '<textFormat foregroundColor="4278190080" IndentMargin="5" Alignment="1" backgroundColor="0" VerticalIndentMargin="50" font="Bahnschrift SemiBold,16,-1,5,50,0,0,0,0,0"/>'
GTextAreaTemplate_str += "</text> </area>"


def modification_date(filename_str):
    t = os.path.getmtime(filename_str)
    return datetime.datetime.fromtimestamp(t)


walkThroughFolder_b = False
imageFiles_filename = "imageFiles.txt"
desiredYear = 2022

imageFiles_listPerMonth = {}
for m in range(0, 12):
    imageFiles_listPerMonth[m + 1] = list()

if walkThroughFolder_b:

    startFolder_str = "C:/Users/Wunstorf/Eigene Bilder/**"
    lastPrintLength = 0
    folderFound_b = False
    for filename_str in glob.iglob(startFolder_str, recursive=True):
        if "jif" in filename_str:
            continue
        if filename_str.endswith("MOV") or filename_str.endswith("mov"):
            continue

        if os.path.isfile(filename_str):  # filter dirs
            date_dtime = modification_date(filename_str)
            if date_dtime.year == desiredYear:
                imageFiles_listPerMonth[date_dtime.month].append(filename_str)
                folderFound_b = True
                # if len(imageFiles) % 20 == 0:
                #      print(f"{len(imageFiles)} : {filename_str}")
        else:
            if folderFound_b:
                print("  --")
            else:
                back_str = "\b" * lastPrintLength
                print(back_str, end="")
            folderFound_b = False

            print(filename_str, end="")
            lastPrintLength = len(filename_str)

    print(" ")

    imageFile_listStr = [
        f"{key},{imageStr}\n"
        for key, imageList in imageFiles_listPerMonth.items()
        for imageStr in imageList
    ]
    with open(imageFiles_filename) as imageFile:
        imageFile.writelines(imageFile_listStr)
else:

    with open(imageFiles_filename) as imageFile:
        imageFile_listStr = imageFile.readlines()

    for csvImageFile_str in imageFile_listStr:
        month_imageFile_list = csvImageFile_str.split(",")
        monthFromStr = int(month_imageFile_list[0])
        imageFilename = month_imageFile_list[1]
        if imageFilename[-1] == "\n":
            imageFilename = imageFilename[0:-1]
        imageFiles_listPerMonth[monthFromStr].append(imageFilename)

for m in range(0, 12):
    print(f"Month:{m} -- imageCnt:{len(imageFiles_listPerMonth[m+1])}")

templateLeerFile_str = "TemplateLeer_3x2.mcf"
with open(templateLeerFile_str) as mcfFile:
    mcfLines_listStr = mcfFile.readlines()

widthHorizontal_str = f'width="{GWidthHorizontal}"'
heightHorizontal_str = f'height="{GHeightHorizontal}"'
imageIdx_list = [
    lineIdx
    for lineIdx, xmlStr in enumerate(mcfLines_listStr)
    if widthHorizontal_str in xmlStr and heightHorizontal_str in xmlStr
]
imageCount = len(imageIdx_list)

destinationJahrbuch_str = f"jahrgangsbuch{desiredYear}"

# Bild laden mit PIL, größe skalieren auf breite, kopieren in dest-folder
# mcf als xml mit template überschreiben
GDestinationFolder_path = Path(f"C:/Users/Wunstorf/Eigenes Photobuch")
GDestinationImageFolder_path = (
    GDestinationFolder_path / f"{destinationJahrbuch_str}_mcf-Dateien"
)
try:
    GDestinationImageFolder_path.mkdir()
except FileExistsError:
    print("Folder is already there")

GDestinationMCF_path = GDestinationFolder_path / f"{destinationJahrbuch_str}.mcf"

selectedImageFilenames_list = []

allImageFiles = [
    filename
    for monthKey, imagefilename_list in imageFiles_listPerMonth.items()
    for filename in imagefilename_list
]
for imageIndex in range(0, imageCount):
    selectedFilename_str = allImageFiles.pop(random.randint(0, len(allImageFiles) - 1))
    selectedImageFilenames_list.append(selectedFilename_str)


random.shuffle(selectedImageFilenames_list)
for imageIdx, selectedFilename_str in enumerate(selectedImageFilenames_list):
    selectedImage_img = Image.open(selectedFilename_str)
    selectedImage_width = selectedImage_img.size[0]
    selectedImage_height = selectedImage_img.size[1]
    scale_width = GWidthHorizontal / selectedImage_width
    scale_height = GHeightHorizontal / selectedImage_height

    selectedFilename_output_path = (
        GDestinationImageFolder_path / "img_{0:06d}.jpeg".format(imageIdx)
    )
    selectedImage_img.save(selectedFilename_output_path, format="jpeg")

print("hut")
