3D-dots-in-elliptical-volume
===

3D dots simulation in a specified elliptical volume.

Use **symula2d.m** to uniformly distribute 2D dots in a given volume.  
Use **symula3d.m** to uniformly distribute 3D dots in a given volume.

Both symula2d.m and symula3d.m are classes. Use `doc symula` and `doc symula3d` for more details and to get an example of how to run the code.

---

**symulaP.lib.R** is a library, it contains the functions to uniformly distribute Cartesian points in a specified elliptical volume.

**symula2d.script.R** is a script, it contains the functions to uniformly distribute 3D dots in a specified elliptical volume. Change the parameters in the top-part of the file and then run it.

## Requirments

### For the R script

To run the **symula2d.script.R** code please make sure to have the `parallel` and `tiff` packages installed in your R environment. You can check that easily by running `library(package_name)` in R. If an error occurs saying that the package cannot be located or is not installed, please install it by running `install.packages("package_name")`.