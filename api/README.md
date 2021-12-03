## Cahaba: Flood Inundation Mapping for U.S. National Water Model

Flood inundation mapping software is configured to work with the U.S. National Water Model operated and maintained by the National Oceanic and Atmospheric Administration (NOAA) National Water Center (NWC).

This software uses the Height Above Nearest Drainage (HAND) method to generate Relative Elevation Models (REMs), Synthetic Rating Curves (SRCs), and catchment grids. This repository also includes functionality to generate flood inundation maps (FIMs) and evaluate FIM accuracy.

## API Tool

This tool allows for batch processing of HUCS.


(More details to come)

### Description:

- This system has three primary components.
	1) FrontEnd
		- The `Gui` (web front end) which is a docker container and another container called the `output_handler`. The Output Handler sends message from the GUI to other containers for processing.
	2) Node
		- This components also has two primary docker containers. One is called `connector` which is for passing data and messages to other locations, folders, or docker containers
		such as the FrontEnd Output Handler.
	3) Nginx
		- Which serves as a webservice as well as use socket.io to help communication between docker containers.

The three components do not necessarily  need to be all running on one host machine. It can be configured to have the "FrontEnd" docker containers running on one host server and the "Node" docker containers running on a seperate host server. This can allow for better resource management and performance if required.
		
### Notes:
	- Some of the .conf, and .yml files have the phrase `-dev` or `-prod` in the file name. When the file name has `-dev`, it means the three components were being deployed all on one host machine. When the file name has `-prod` in it, it means the entire system will be running against two host servers. The "FrontEnd" docker containers will have docker containers for the gui, output hander and an nginx image.  The "Node" docker containers will have docker containers for connector, updater and its own nginx image.