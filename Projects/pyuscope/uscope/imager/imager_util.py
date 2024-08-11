from PIL import Image
try:
    import usb
except ImportError:
    usb = None

if usb:
    import usb.core
    import usb.util
import glob
import os
import subprocess
from uscope.util import tostr

VID_TOUPTEK = 0x0547
PID_MU800 = 0x6801


def get_scaled(image, factor, filt=Image.NEAREST):
    if factor == 1.0:
        return image
    else:
        return image.resize(
            (int(image.size[0] * factor), int(image.size[1] * factor)), filt)


def have_touptek_camera():
    if usb:
        # find our device
        for dev in usb.core.find(find_all=True):
            if dev.idVendor != VID_TOUPTEK:
                continue
            return True


def have_v4l2_camera():
    # FIXME: more proper check
    return os.path.exists("/dev/video0")


def auto_detect_source(verbose=False):
    # FIXME: this is probably obsolete now
    # Way to detect if it is a modifed driver w/ direct gain control?
    # Takes about 20 ms, leave in
    if b"touptek" in subprocess.check_output("lsmod"):
        for dev in usb.core.find(find_all=True):
            if dev.idVendor != VID_TOUPTEK:
                continue
            verbose and print(
                "ADS: found ToupTek kernel module + camera (MU800?)")
            assert glob.glob("/dev/video*"), "Camera not found???"
            return "gst-v4l2src-mu800"

    if usb:
        for dev in usb.core.find(find_all=True):
            if dev.idVendor != VID_TOUPTEK:
                continue
            verbose and print("ADS: found VID 0x%04X, PID 0x%04X" %
                              (dev.idVendor, dev.idProduct))
            verbose and print("ADS: found ToupTek generic camera")
            return "gst-toupcamsrc"

    # Fallback to generic gst-v4l2 if there is an unknown camera
    if glob.glob("/dev/video*"):
        verbose and print("ADS: found /dev/video device")
        return "gst-v4l2src"

    verbose and print("ADS: giving up (usb: %s)" % bool(usb))
    return "gst-testsrc"


def format_mm_3dec(value):
    """
    Given a value in mm return a string printing it nicely
    Currently always gives 3 decimal places
    """
    if value >= 100:
        return "%d mm" % value
    elif value >= 10:
        return "%0.1f mm" % value
    elif value >= 1:
        return "%0.2f mm" % value
    elif value >= 0.1:
        return "%d um" % (value * 1000, )
    elif value >= 0.01:
        return "%0.1f um" % (value * 1000, )
    elif value >= 0.001:
        return "%0.2f um" % (value * 1000, )
    else:
        if value < 0.0001 and value != 0.0:
            print("WARNING: formatitng very small number: %g" % (value, ))
        return "%0.3f um" % (value * 1000, )