Phantom Viewer
==============

A short tutorial.

Cameras
-------

- `phantom2`: divertor
- `phantom`: outboard midplane

Play videos
-----------

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

Get video data
--------------

    import acquire
    time = acquire.gpi_series(shot, 'phantom2', 'time')
    frames = acquire.video(shot, camera='phantom2', sub=20, cache=False)

- `cache=True` saves a numpy object of the video to the `cache` directory

Reconstruct a video
-------------------

    from phantom_viewer.fl import reconstruct
    reconstruct.write_nnls_reconstruction(1150625010)

Takes about an hour.

View a reconstruction
---------------------

    from phantom_viewer.fl import view
    view.slide_reconstruction(1150625010)

Setup
-----

Make sure you run the code on a Fedora machine at the PSFC, the Red Hat ones have a different Python installation that has problems.

### Just want to browse data

    cd /home/sballinger/Desktop/phantom_viewer
    source venv/bin/activate

### Need to get data for your own scripts

In the directory with your scripts, do

    ln -s /home/sballinger/Desktop/phantom_viewer phantom_viewer

then import all components in your script like this:

    from phantom_viewer import blah

### Make your own copy of the code

    git clone https://github.com/sballin/phantom_viewer
    cd phantom_viewer

Write scripts here or run the above commands from the console. If it complains about missing packages, do `pip install --user [package name]`.
