#!/usr/bin/env python3
"""
Running the full suite:
-GRBL controller attached (no microscope)
-E3ISPM20000KPA camera attached
-v4levice as /dev/video0 that supports 640x480 video
    Ex: my X1 carbon has this as built in web camera
"""

import unittest
import os
import json5
from uscope.motion.hal import MockHal
from uscope.imager.imager import MockImager
from uscope.imager import gst
from uscope.imager.util import have_touptek_camera, have_v4l2_camera
import uscope.planner
from uscope.config import get_usj
from uscope.util import printj
from uscope.planner import microscope_to_planner
import shutil
import time
import glob
from uscope.motion.grbl import GRBL, GrblHal
from uscope.imager.touptek import toupcamsrc_info

TT_WH = (5440, 3648)
TT_WH = (4928, 4928)
V4L2_WH = (640, 480)


class TestCommon(unittest.TestCase):
    def setUp(self):
        """Call before every test case."""
        print("")
        print("")
        print("")
        print("Start " + self._testMethodName)
        self.verbose = os.getenv("VERBOSE", "N") == "Y"
        self.verbose = int(os.getenv("TEST_VERBOSE", "0"))
        self.planner_dir = "/tmp/pyuscope/planner"
        if os.path.exists("/tmp/pyuscope"):
            shutil.rmtree("/tmp/pyuscope")
        os.mkdir("/tmp/pyuscope")

    def tearDown(self):
        """Call after every test case."""


class PlannerTestCase(TestCommon):
    def simple_planner(self, pconfig, motion=None, imager=None, dry=False):
        if self.verbose:
            log = print
        else:
            log = lambda x: None
        if motion is None:
            motion = MockHal(log=log)
        if imager is None:
            imager = MockImager(width=150, height=50)
        planner = uscope.planner.Planner(pconfig,
                                         motion=motion,
                                         imager=imager,
                                         out_dir=self.planner_dir,
                                         progress_cb=None,
                                         dry=dry,
                                         log=log,
                                         verbosity=2)
        return planner.run()

    def simple_config(self):
        """
        Simple scan config
        4 images wide
        3 images tall
        """
        return {
            # "microscope": json5.load(open("configs/mock/microscope.j5", "r")),
            # "motion": {
            #    "origin": "ll",
            #},
            "imager": {
                #    "scalar": 0.5,
                "x_view": 1.0,
            },
            "contour": {
                "start": {
                    "x": 0.0,
                    "y": 0.0,
                },
                "end": {
                    "x": 2.0,
                    "y": 1.0,
                },
            },
        }

    def test_dry(self):
        """blah"""
        meta = self.simple_planner(pconfig=self.simple_config(), dry=True)
        # printj(meta)
        expect_images = 12
        self.assertEqual(expect_images, meta["planner"]["pictures_taken"])
        self.assertEqual(expect_images, meta["planner"]["pictures_to_take"])
        self.assertEqual(expect_images, len(meta["images"]))

    def test_simple(self):
        meta = self.simple_planner(pconfig=self.simple_config(), dry=False)
        # printj(meta)
        expect_images = 12
        self.assertEqual(expect_images, meta["planner"]["pictures_taken"])
        self.assertEqual(expect_images, meta["planner"]["pictures_to_take"])
        self.assertEqual(expect_images, len(meta["images"]))

    def test_tsettle(self):
        """
        test pconfig["tsettle"]
        This controls how long to wait between movement and snapping a picture
        """

        pconfig = self.simple_config()
        pconfig["tsettle"] = 0.0
        # get a baseline without
        tstart = time.time()
        self.simple_planner(pconfig=pconfig)
        d0 = time.time() - tstart

        pconfig["tsettle"] = 0.01
        tstart = time.time()
        meta = self.simple_planner(pconfig=pconfig)
        d1 = time.time() - tstart

        # should be 12 images and added 0.1 sec per image
        # so should have increased by at least a second
        # make it quicker...I'm impatient
        self.assertEqual(12, meta["planner"]["pictures_taken"])
        d = d1 - d0
        self.verbose and print("delta: %0.2f" % d)
        assert d > 0.1

    def test_scalar(self):
        pconfig = self.simple_config()
        pconfig["imager"]["scalar"] = 0.5
        self.simple_planner(pconfig=pconfig)

    def test_backlash(self):
        pconfig = self.simple_config()
        pconfig.setdefault("motion", {})["backlash"] = 0.1
        self.simple_planner(pconfig=pconfig)

    def test_origin_ll(self):
        pconfig = self.simple_config()
        pconfig.setdefault("motion", {})["origin"] = "ll"
        self.simple_planner(pconfig=pconfig)

    def test_origin_ul(self):
        pconfig = self.simple_config()
        pconfig.setdefault("motion", {})["origin"] = "ul"
        self.simple_planner(pconfig=pconfig)

    def test_coordinates_rounded(self):
        """
        See if pictures are taken at the right coordinates
        There should be some tolerance when very cloose
        """
        pconfig = {
            "step": 0.75,
            "imager": {
                "x_view": 1.0,
                # y_view => 0.5
            },
            "contour": {
                "start": {
                    "x": 0.0,
                    "y": 0.0,
                },
                "end": {
                    # Exactly 3 images wide
                    # 1.0 * (1 + 2 * 0.75)
                    "x": 2.5,
                    # Exactly 2 images tall
                    # 0.5 * (1 + 1 * 0.75)
                    "y": 0.875,
                },
            },
        }
        imager = MockImager(width=100, height=50)
        j = self.simple_planner(imager=imager, pconfig=pconfig)
        # y = 0.5 * 0.75 = 0.375
        expect_coordinates = [
            (0, 0),
            (0.75, 0),
            (1.5, 0),
            (1.5, 0.375),
            (0.75, 0.375),
            (0, 0.375),
        ]
        self.validate_coordinates(j, expect_coordinates)

    def validate_coordinates(self, j, expect_coordinates):
        """
        Dict order in order taken
        Can be sorted to row/column order
        """
        self.assertEqual(len(expect_coordinates), len(j["images"]))
        # JSON outputs each column at a time regardless of order taken
        for gotj, expect in zip(j["images"].values(), expect_coordinates):

            def xy_delta(a, b):
                return ((a[0] - b[0])**2 + (a[1] - b[1])**2)**0.5

            delta = xy_delta(expect, (gotj["x"], gotj["y"]))
            self.assertTrue(delta < 0.01, (expect, gotj, delta))

    def test_exclude(self):
        pconfig = self.simple_config()
        pconfig["exclude"] = [{"r0": 0, "c0": 0, "r1": 1, "c1": 1}]
        self.simple_planner(pconfig=pconfig)

    def test_microscope_to_planner(self):
        usj = get_usj(name="mock")
        contour = {
            "start": {
                "x": 0.0,
                "y": 0.0,
            },
            "end": {
                "x": 2.0,
                "y": 1.0,
            }
        }
        microscope_to_planner(usj, objectivei=0, contour=contour)


class GstTestCase(TestCommon):
    def test_mock(self):
        usj = get_usj(name="mock")
        gst.get_cli_imager_by_config(usj)


class GstCLIImagerTestCase(TestCommon):
    def get_image(self):
        ret = []

        def thread(loop):
            # self.imager.warm_up()
            im = self.imager.get()
            ret.append(im)

        self.imager.gst_run(thread)
        if len(ret) == 0:
            raise Exception("No images")
        return ret[0]

    def test_get_args(self):
        import argparse

        parser = argparse.ArgumentParser(
            description="GstImager (gstreamer wrapper) demo")
        gst.gst_add_args(parser)
        args = parser.parse_args([])

        gst.gst_get_args(args)

    def test_gst_failed_start(self):
        """
        If gstreamer gets a bad device it should shut down gracefully
        """
        # FIXME: only run if there is an appropriate device
        self.imager = gst.GstCLIImager(
            opts={
                "source": "v4l2src",
                "v4l2src": {
                    "device": "/dev/video123",
                },
                "wh": (640, 480)
            })
        try:
            _images = self.get_image()
            self.fail("Expected exception")
        except gst.GstError:
            pass

    def test_videotestsrc(self):
        self.imager = gst.GstCLIImager(opts={
            "source": "videotestsrc",
            "wh": (123, 456)
        })
        im = self.get_image()["0"]
        self.assertEqual((123, 456), im.size)

    def test_raw(self):
        """
        Need 59535360 bytes, got 59535360
        """
        if not have_touptek_camera():
            self.skipTest("no touptek camera")
        # Doesn't work...hmm
        # maybe its jpgeg or something like that?
        # self.imager = gst.GstCLIImager(opts={"source": "videotestsrc", "wh": (100, 100), "gst_jpg": False})
        self.imager = gst.GstCLIImager(opts={
            "source": "toupcamsrc",
        })
        im = self.get_image()["0"]
        self.assertEqual((TT_WH), im.size)

    def test_v4lsrc(self):
        """
        video4linux test
        (using my laptop camera)
        """
        if not have_v4l2_camera():
            self.skipTest("no v4l2 camera")
        # FIXME: only run if there is an appropriate device
        self.imager = gst.GstCLIImager(
            opts={
                "source": "v4l2src",
                "v4l2src": {
                    "device": "/dev/video0",
                },
                "wh": V4L2_WH,
            })
        im = self.get_image()["0"]
        self.assertEqual(V4L2_WH, im.size)

    def test_toupcamsrc(self):
        """
        touptek test
        (using E3ISPM20000KPA (IMX183))
        """
        if not have_touptek_camera():
            self.skipTest("no touptek camera")
        # FIXME: only run if there is an appropriate device
        self.imager = gst.GstCLIImager(opts={
            "source": "toupcamsrc",
        })
        im = self.get_image()["0"]
        self.assertEqual(TT_WH, im.size)
        self.imager.stop()

    def test_toupcamsrc_info(self):
        if not have_touptek_camera():
            self.skipTest("no touptek camera")
        j = toupcamsrc_info()
        assert j["eSizes"][0]["StillResolution"]["w"] > 0


def get_grbl(verbose=False):
    # print("Checking for GRBL...")
    if not glob.glob("/dev/serial/by-id/usb-*_USB_Serial-if00-port0"):
        print("GRBL: no plausible GRBL")
        return None
    try:
        ret = GRBL(verbose=verbose)
    except:
        print("GRBL: open failed")
        raise
    print("GRBL: open ok")
    return ret


class GrblTestCase(TestCommon):
    def setUp(self):
        """Call before every test case."""
        super().setUp()

        self.grbl = get_grbl()
        if not self.grbl:
            self.skipTest("No GRBL")

        # print("Set TEST_GRBL_ALL to run aggressive tests")

        # A small amount to move
        self.delta = 0.01
        # A reasonable feedrate
        self.f = 100

    def tearDown(self):
        """Call after every test case."""
        self.grbl.close()
        super().tearDown()

    def test_reset(self):
        self.grbl.reset()

    def test_stop(self):
        self.grbl.stop()

    def test_qstatus(self):
        self.grbl.qstatus()

    def test_mpos(self):
        self.grbl.mpos()

    def test_move_absolute(self):
        mpos = self.grbl.mpos()
        for axis in "xyz":
            pos0 = mpos[axis]
            self.grbl.move_absolute(pos={axis: pos0 + self.delta}, f=self.f)
            self.grbl.move_absolute(pos={axis: pos0}, f=self.f)

    def test_hard_move_relative(self):
        for axis in "xyz":
            self.grbl.move_relative(pos={axis: +self.delta},
                                    f=self.f,
                                    soft=False)
            self.grbl.move_relative(pos={axis: -self.delta},
                                    f=self.f,
                                    soft=False)

    def test_soft_move_relative(self):
        for axis in "xyz":
            self.grbl.move_relative(pos={axis: +self.delta},
                                    f=self.f,
                                    soft=True)
            self.grbl.move_relative(pos={axis: -self.delta},
                                    f=self.f,
                                    soft=True)

    def test_wait_idle(self):
        self.grbl.wait_idle()

    def test_jog(self):
        self.grbl.jog({"x": +self.delta}, rate=self.f)
        time.sleep(0.1)
        self.grbl.cancel_jog()
        self.grbl.jog({"x": -self.delta}, rate=self.f)
        time.sleep(0.1)
        self.grbl.cancel_jog()

    def test_cancel_jog(self):
        self.grbl.cancel_jog()


class GrblHalTestCase(TestCommon):
    def setUp(self):
        """Call before every test case."""
        super().setUp()

        grbl = get_grbl()
        if not grbl:
            self.skipTest("No GRBL")
        self.gh = GrblHal(grbl=grbl)

        # print("Set TEST_GRBL_ALL to run aggressive tests")

        # A small amount to move
        self.delta = 0.01
        # A reasonable feedrate
        self.f = 100

    def tearDown(self):
        """Call after every test case."""
        self.gh.close()
        super().tearDown()

    def test_pos(self):
        self.gh.pos()

    def test_move_absolute(self):
        pos0s = self.gh.pos()
        for axis in "xyz":
            pos0 = pos0s[axis]
            self.gh.move_absolute(pos={axis: pos0 + self.delta})
            self.gh.move_absolute(pos={axis: pos0})

    def test_move_relative(self):
        for axis in "xyz":
            self.gh.move_relative(pos={axis: +self.delta})
            self.gh.move_relative(pos={axis: -self.delta})


if __name__ == "__main__":
    unittest.main()
