SIMPLE script for batch mesh rendering.

FILES:
	- blenderBatchRender.py : 
			Python script for executing the batch processing. 
			Takes as input 3 arguments: <inputDIR> <globPattern> <outputDIR>
			<inputDIR> contains the input meshes
			<globPattern> glob pattern used to retrieve the meshes inside <inputDIR>
			<outputDIR> directory to store the snapshot images
					
	- snapshots.blend : 
			Blender project file used to setup scene related stuff, like lighting, materials, base scene etc..

USAGE:

	- Setup your scene tuning snapshots.blend, then call the script like this:
		blender --background snapshots.blend -P blenderBatchRender.py -- <directoryIN> <globPattern> <directoryOUT>