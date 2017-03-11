#####<a href="../README.html">&#9664; Projects index</a>

ND2 to TIFF file conversion
===

## Convert_to_tif.js

`Convert_to_tif.js` is a JavaScript script that allows to convert images from the `.nd2` (and other formats) to the `.tif` format. To use it:

1. Open FIJI
2. Drag and drop the `Convert_to_tif.js` file on FIJI
3. Click on the `Run` button.
4. Select one of the `.nd2` files in the folder and press OK. This step is required just to identify the folder. Please note that **EVERY** `.nd2` file in that folder will be converted.
    - The script detects the extension of the files you want to convert from the selected one, and then convert all the file sin that folder with that extension. For example, if you select a `.czi` file in a folder, the script will convert all `czi` files in that folder to `.tif`.
5. Type the channel names separated by commas (no white-spaces) and press OK. Please note that the channels <u>must be in the proper order</u>.
6. Select the file-name notation to be used for the output:
    - DOTTER: dapi_001.tif
    - GPSeq: dapi.channel001.series001.tif
7. The script will create a subfolder for every nd2 file in the selected folder and will create single-channel tif images with the chosen name notation.

An updated FIJI is required to run the script.

#### How to permanently add the script to Fiji

Copy the `Convert_to_tif.js` file to your `./Fiji.app/plugins/Scripts/Image/` folder and rename it as 'Convert_to_tif.js'. If you restart FIJI you'll see a new option in the `Image` menu.

## nd2toTif.py

**nd2toTif.py** is a Python script that allows to convert images from the `.nd2` to the `.tif` format. The script can accept either the folder where the nd2 files are located or the path to the nd2 file.

    ./nd2toTif.py /home/usr/nd2-file-folder/
    ./nd2toTif.py /hom/user/file.nd2

#### Requirements

Please make sure to have the `javabridge`, `bioformats` and `numpy` packages installed in your Python environment. You can check that easily by running `import package_name` in Python. If an error occurs saying that the package cannot be located or is not installed, please install it by running `sudo pip install package-name` or `sudo python -m pip install package-name`.
