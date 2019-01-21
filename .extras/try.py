import Imagetree as it
import cv2
import numpy as np 
from imutils.video import VideoStream
import imutils
import matplotlib.pyplot as plt
import matlplotlib.image as mpimg
import time


#Defining the color of the object in the frame
#Color  = Green  ====> Change to red
colorLower  = it.Pixel(29,86,6)
colorUpper = it.Pixel(64,255,255)

#Starting the videostream

vs  = VideoStream(src=0).start()

time.sleep(2.0)

while True:
    frame =  vs.read()
    if frame is None:
        break
    frame = frame
    frame = imutils.resize(frame,width=600)
    key = cv2.waitKey(1)
    blurred = cv2.GaussianBlur(frame,(11,11),0)
    mask = cv2.inRange(hsv,colorcolorLower)
    mask = cv2.erode(mask,None,iterations=2)
    mask = cv2.dilate(mask,None,iterations=2)
    cv2.imshow("Frame",mask)
    key = cv2.waitkey(1)



















# leaf1 =  it.QuadLeaf()
# print(leaf1)
# pix1 = it.Pixel(0,0,255,0)
# leaf1 = it.QuadLeaf(pix1)
# print(leaf1)
# print(pix1)
# pix1 = it.Pixel(255,0,255,0)

