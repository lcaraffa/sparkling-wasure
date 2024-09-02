# 1.7.2
- Add possibility to select extra dimensions to keep in standardization

# 1.7.1
Same as 1.7.0 (new tag needed to publish on pypi due to incorrect package handling)

# 1.7.0
- las_remove_dimension: new tool to remove one or many dimensions
- deploy on ghcr.io instead of dockerhub
- Add tools to run functions on buffered las:
  - update create_las_with_buffer to enable saving which points are from the central las on a new dimension
  - add a remove_points_from_buffer to remove the points that have this new dimension not set to 1
  - add a decorator to run a function on a buffered las and return an output las only with the points from the original input

# 1.6.0
- color: choose streams for RGB colorization, and IRC colorization (doc https://geoservices.ign.fr/services-web-experts-ortho)
- color: detect white images.

# 1.5.2
- refactor tool to propagate header infos from one pipeline to another to use it by itself

# 1.5.1
- fix add_buffer: propagate header infos from input to the output
- update pdal.Writer params to make sure input format is forwarded except for the specified parameters
- add test for colorization with epsg != 2154

# 1.5.0
- switch colorisation source from Geoportail to Geoplateforme
- use absolute value comparison in tests
- linting / add pre-commits
- upgrade python version to 3.11
- replace `pdal info --metadata` subprocess by a python command (used in the `color` module in particular)

# 1.4.1
- fix copy_and_hack_decorator (was not returning the decorated function output)

# 1.4.0
- count_occurences / replace_value: add copy_and_hack decorator to run on tscan output files
- Update to pdal 2.6+ to better handle classification values and flags in replace_attribute_in_las
(was treating values over 31 as {classification under 31 + flag} even when saving to LAS 1.4)

# 1.3.1
- fix color: ensure that tmp orthoimages are deleted after use by using the namedTemporaryFile properly.

# 1.3.0
- color: support colorization for <0.2m clouds (including height=0/width=0)
- color: ceil width/height to have a bbox that contains all points

# 1.2.1
- fix cicd_full github action: deployment was triggered on pushing to dev instead of master only

# 1.2.0
- color: keep downloaded orthoimages by returning them to make them stay in execution scope

# 1.1.1
- unlock: fix main
- tests:
  - add geoportail marker to skip tests relying on geoportail (they are now played on PR to the master branch only)
  - bugfix on standardization test file path

# 1.1.0
- standardization: handle malformed laz input ("Global encoding WKT flag not set for point format 6 - 10")
- color: extract unlock module from colorization and rename colorization function

# 1.0.0
- first public version
- docker: Use staged build to reduce docker image size
- add continuous integration with github actions

# v0.5.6
- makefile: run "clean" before "build" (build the library) in order for the CI to remove old versions of the library

# v0.5.5
- standardisation: set offset to 0

# v0.5.4
- standardisation: fix warnings displayed when using lasinfo (LasTools). Use las2las to save the las again

# v0.5.3
- add_buffer/merge: use filename to get tile extent

# v0.5.2 :
- jenkins script: handle errors
- docker: inherit of an image based on Mamba instead of Conda (mamba is faster than conda to fetch dependencies)
- continuous integration (jenkins): build the docker image and publish on our private repo when merging on master

# v0.5.1
- standardisation : parallelize occurences count

# v0.5.0
- docker: option no-capture-output
- standardisation : add a module to enforce format for a las/laz file
- standardisation : add a module to count occurences for an attribute in a batch of las/laz file
- standardisation : add a module to replace the values of an attribute in a las/laz file

# v0.4.2
standardisation
stitching
