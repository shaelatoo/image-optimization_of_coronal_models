# image-optimization_of_coronal_models
This code is used to create image-optimized models of the solar coronal magnetic field, based on coronagraph images.

The model is constructed using the PFSS SolarSoft library created by Marc DeRosa.  This code in the library 
extrapolates the coronal magnetic field from an input synoptic magnetogram.  An unfinished but fairly thorough explanation of the
process, the programs included in the library, and the variables used by the program can be found at:
http://www.lmsal.com/~derosa/pfsspack/ and multiple explanations of the mathematics of how this is traditionally achieved can be
found in the literature, (e.g. https://link.springer.com/content/pdf/10.12942%2Flrsp-2012-6.pdf), though see
http://iopscience.iop.org/article/10.1088/0004-637X/732/2/102/meta for an excellent discussion of some problems with the traditional
approach that may make a numerical solution preferable in some situations.

The primary function of the optimization code is to compare data derived from coronagraph observations to the model, calculate
an objective function that quantifies the level of disagreement between the two, and then iteratively alter the spherical 
harmonic coefficients of the original synoptic magnetogram via the Downhill Simplex optimization algorithm to achieve better
agreement between the model and the image-based constraints.  

What you might call "Main" is the function harmonic_amoeba2.pro, which contains the optimization algorithm.  The repository also
contains the following programs:

pfss_opt_parameters.pro - contains a common block variable that contains parameters used by the PFSS library and the optimization
    software
transform_coefficients.pro - contains a common block variable that contains stored values related to the re-calculation of
    the inverse harmonic transforms by optimization_inv_transform.pro
calc_phis.pro - given the modified values of the magnetogram transform, updates the values of phiat and phibt in the PFSS common
    block
extract_transform.pro - given a vector of spherical harmonic transform coefficient values (a vertex in the optimization scheme),
    reconstructs the 2D harmonic transform
find_angle_values2 - samples the magnetic field model at specified locations and projects it onto the corresponding image plane,
    returning the orientation angle (0 to pi) of the projected field in the plane
form_initial_vertex.pro - forms the initial vertex from the original harmonic transform; only certain harmonic coefficients are 
    optimized, so these are extracted from the transform and put into a vector; this routine and extract_transform.pro are 
    essential inverse of one another
harmonic_trypoint2.pro - calculates the value of the objective function for a given vertex; per amoeba
    implementation of the downhill simplex method, it may alter the simplex, but it can be called for penalty calculation only,
    per the penalty_only keyword
inv_spherical_transform_sj.pro - an altered version of Marc DeRosa's inv_spherical_transform.pro; this altered version was 
    created to take advantage of the fact that only part of the inverse transforms need to be recalculated during each iteration 
    of the optimization algorithm
optimization_inv_transform.pro - this is the routine that re-calculates the inverse transforms once they've been calculated; it 
    interacts with inv_spherical_transform_sj.pro
pfss_potl_field_sj.pro - modified version of Marc DeRosa's pfss_potl_field.pro; modification created to use my versions of the 
    inverse spherical harmonic transform; may contain some other small updates
pfss_mag_create_sj.pro - my version of Marc DeRosa's pfss_mag_create.pro; initially created to correct an issue I was having with
    hourly GONG magnetograms (they needed to be shifted so that the leftmost pixel column corresponded to carrington longitude
    zero), this routine now includes code that remove net flux from the magntograms before interpolating on the PFSS model grid
    and code that allows them to work with ADAPT ensemble maps
zero_array.pro - performs net flux removal for pfss_mag_create_sj.pro
harmonic_component_variability.pro - calculates the range of values for each spherical harmonic transform coefficient over a set
    of synoptic magnetogram maps; used to set scale variable for harmonic_amoeba2.pro
find_magtype.pro - magtype is a variable used by PFSS software library to determine how to read the synoptic magnetograms; this
    routine should determine it, but currently only discriminates between magfile=2 and magfile=4
