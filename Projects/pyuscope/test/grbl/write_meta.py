#!/usr/bin/env python3
from uscope.motion.grbl import GRBL, grbl_write_meta, grbl_read_meta, NoGRBLMeta, microscope_name_hash, grbl_delete_meta
from uscope.util import add_bool_arg
import os


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Write pyuscope metadata to GRBL controller")
    add_bool_arg(parser, "--verbose", default=False, help="Verbose output")
    parser.add_argument("--comment", help="Comment. 8 chars max")
    parser.add_argument("--sn", help="Serial number. 8 chars max")
    parser.add_argument("--microscope",
                        help="Microscope config file name to hash")
    add_bool_arg(parser, "--delete", default=False, help="Delete metadata")
    args = parser.parse_args()

    grbl = GRBL()

    if args.delete:
        print("Deleting")
        grbl_delete_meta(grbl.gs)
    else:
        comment = args.comment
        sn = args.sn

        microscope = args.microscope
        config = None
        if microscope:
            microscope_dir = os.path.join("configs", microscope)
            if not os.path.exists(microscope_dir):
                raise Exception("Invalid microscope")
            # sha256 => resulted in error
            # maybe not quite enough bits?
            # just roll with this and move on
            # config = hashlib.sha256(microscope.encode("ascii")).digest()[0:4]
            config = microscope_name_hash(microscope)
            print("Microscope %s => config %s" % (microscope, config.hex()))

        info = {}
        if comment is None or sn is None or config is None:
            try:
                info = grbl_read_meta(grbl.gs)
                print("Queried existing info")
                print("  Comment: %s" % (info["comment"], ))
                print("  S/N: %s" % (info["sn"], ))
                print("  Config: %s" % (info["config"].hex(), ))
            except NoGRBLMeta:
                print("WCS magic number not found")
        if sn is None:
            sn = info.get("sn")
        if comment is None:
            comment = info.get("comment")
        if config is None:
            config = info.get("config")
        print("Writing")
        grbl_write_meta(grbl.gs, comment=comment, sn=sn, config=config)

    print("Verifying")
    try:
        info = grbl_read_meta(grbl.gs)
    except NoGRBLMeta:
        print("Config magic number not found")
        return
    print("  Comment: %s" % (info["comment"], ))
    print("  S/N: %s" % (info["sn"], ))
    print("  Config: %s" % (info["config"].hex(), ))
    assert info["comment"] == comment
    assert info["sn"] == sn
    assert info["config"] == config

    print("Write ok")


if __name__ == "__main__":
    main()
