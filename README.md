Phantom Viewer
==============

A short tutorial.

Cameras
-------

- `phantom2`: divertor
- `phantom`: outboard midplane

Playing videos
--------------

    import view

Explore with a slider:

    view.slide_phantom(1150625010) # can take up to 30 seconds

- Click on the `Sub 20` option (and wait for a while) to do background subtraction
- Click the `Recolor` button to set the min/max colors to the min/max values of the currently displayed image
- Use the slider to move around the video
- Use the `>` and `<` buttons to move one frame forward/backward
- Click play to play 100 frames

Just play the whole thing at a good speed:

    view.animate_phantom(1150625010, sub=20, interval=10, skip=1)

- All named arguments in this tutorial are defaults and can be omitted
- `sub=0` gets the original video
- `skip=N` will display every Nth frame
- `interval=20` and greater will play the video more slowly

Getting video data
------------------

    import acquire
    time = acquire.gpi_series(shot, 'phantom2', 'time')
    frames = acquire.video(shot, camera='phantom2', sub=20, cache=False)

- `cache=True` saves a numpy object of the video to the `cache` directory

Installing
----------

Clone this repo, `cd` inside, and write scripts there or run the above commands from the console. If it complains about lacking packages, do `pip install --user [package name]` on PSFC machines.

