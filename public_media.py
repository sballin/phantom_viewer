import os
import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']='PySide'
import matplotlib.pyplot as plt
import numpy as np
import acquire
import process


def fit_frame():
    shot = 1150717011
    frames = acquire.video(shot, 'phantom2', sub=20)
    fig = plt.figure(frameon=False)
    fig.set_size_inches(1,1)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    i = 1
    for frame in [6172]:
        plt.imshow(frames[frame], cmap=plt.cm.gray, origin='bottom', aspect='normal')
        plt.autoscale()
        fig.savefig('out/fit_%d.pdf' % i)
        i += 1
    

def elm_frames():
    shot = 1150805004
    frames = acquire.video(shot, 'phantom2', sub=20)
    fig = plt.figure(frameon=False)
    fig.set_size_inches(1,1)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    i = 1
    for frame in [38214, 38218, 38222, 38226, 38230]:
        plt.imshow(frames[frame], cmap=plt.cm.hot, origin='bottom', aspect='normal')
        plt.autoscale()
        fig.savefig('out/elm_%d.pdf' % i)
        i += 1


def single_double_frames():
    shot = 1150730024
    frames = acquire.video(shot, 'phantom2')
    frames_sub = acquire.video(shot, 'phantom2', sub=20)
    fig = plt.figure()
    i = 1
    for frame in [5769, 8413]:
        plt.imshow(frames[frame], cmap=plt.cm.hot, origin='bottom', aspect='normal')
        plt.colorbar()
        plt.autoscale()
        fig.savefig('out/single_double_nosub_%d.pdf' % i)
        plt.clf()
        i += 1


def LH_specgram():
    shot = 1150625025
    frames = acquire.video(shot, 'phantom2', sub=20)
    time_camera = acquire.gpi_series(shot, 'phantom2', 'time')
    time_step = (time_camera[-1]-time_camera[0])/float(frames.shape[0])
    time_dens, dens = acquire.time_dens(shot)
    nfft = 512
    fig, ax = plt.subplots()
    plt.subplot(211)
    plt.specgram(frames[:, 10, 10], NFFT=nfft, window=np.hanning(nfft), detrend='linear', Fs=1./time_step, noverlap=nfft/2)
    plt.xlim([0, time_camera[-1]-time_camera[0]])
    plt.ylabel('Frequency (Hz)')
    plt.subplot(212)
    plt.plot(time_dens, dens)
    plt.xlim([time_camera[0], time_camera[-1]])
    plt.ylim([1.4, 2.5])
    plt.ylabel('Density')
    plt.xlabel('Time (s)')
    plt.show()


def LH_frames():
    shot = 1150625025
    video = acquire.video(shot, 'phantom2', sub=20)
    fig = plt.figure()
    print dir(fig)
    i = 1
    for frame in [4807, 27283]:
        plt.imshow(video[frame], cmap=plt.cm.hot, origin='bottom', aspect='normal')
        plt.colorbar()
        plt.autoscale()
        fig.savefig('out/LH_%d.pdf' % i)
        plt.clf()
        i += 1


def outward_moving_filament_frames():
    shot = 1150717011
    video = acquire.video(shot, 'phantom2', sub=20)
    fig = plt.figure(frameon=False)
    fig.set_size_inches(1,1)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    i = 1
    for frame in [7350, 7354, 7358, 7362]:
        ax.imshow(video[frame], cmap=plt.cm.hot, origin='bottom', aspect='normal')
        plt.autoscale()
        fig.savefig('out/outward_%d.pdf' % i)
        i += 1


def inward_moving_filament_frames():
    shot = 1150618018
    video = acquire.video(shot, 'phantom2', sub=20)
    fig = plt.figure(frameon=False)
    fig.set_size_inches(1,1)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    i = 1
    for frame in [11820, 11824, 11828, 11832, 11836]:
        ax.imshow(video[frame], cmap=plt.cm.hot, origin='bottom', aspect='normal')
        plt.autoscale()
        fig.savefig('out/inward_%d.pdf' % i)
        i += 1


def comet_frames():
    shot = 1150611004
    video = acquire.video(shot, 'phantom2', sub=20)
    fig = plt.figure(frameon=False)
    fig.set_size_inches(1,1)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    i = 1
    for frame in [7450, 7451, 7452, 7453]:
        ax.imshow(video[frame], cmap=plt.cm.hot, origin='bottom', aspect='normal')
        plt.autoscale()
        fig.savefig('out/comet_%d.pdf' % i)
        i += 1


def dens_slices():
    shot = 1150528015
    time, dens = acquire.time_dens(shot)
    plt.figure()
    plt.plot(time, dens)
    plt.axvline(.60074, c='r')
    plt.axvline(.62916, c='r')
    plt.axvline(.66772, c='r')
    plt.axvspan(.6, .67676, color='#dddddd', lw=0)
    plt.xlabel('Time (s)')
    plt.ylabel('Line average density')
    plt.xlim([.5, 1])
    plt.show()


def dens_gif():
    # Also did for frames 5384:6000
    shot = 1150528017
    time, dens = acquire.time_dens(shot)
    time_phantom = acquire.gpi_series(shot, 'phantom2', 'time')
    time, dens = process.time_crop((time, dens), time_phantom)

    plt.figure()
    plt.plot(time, dens)
    plt.xlim([time_phantom[0], time_phantom[-1]])
    plt.xlabel('Time (s)')
    plt.ylabel('Line average density')
    plt.axvspan(time_phantom[28846], time_phantom[29807], color='#dddddd', lw=0)
    for i in range(28846, 29807):
        v = plt.axvline(time_phantom[i], color='r')
        plt.savefig('outframe_%05d.png' % i)
        v.remove()
    os.system('convert outframe_*.png -layers optimize out.gif')
    os.system('rm outframe_*.png')


def output_talk_frames():
    output_frames(1150611004, .761, .762)
    acquire.Database().purge()
    output_frames(1150717011, .752, .7542)
    acquire.Database().purge()
    output_frames(1150625030, .605, .611)
    acquire.Database().purge()
    output_frames(1150820011, 1.2, 1.205)

