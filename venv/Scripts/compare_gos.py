#!C:\Users\Ardopa\Desktop\SBI\Colm_Asia\mi_herramienta_web\venv\Scripts\python.exe
"""Compare two or more sets of GO IDs. Best done using sections."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.cli.compare_gos import CompareGOsCli


def run():
    """Compare two or more sets of GO IDs. Best done using sections."""
    obj = CompareGOsCli()
    obj.write(obj.kws.get('xlsx'), obj.kws.get('ofile'), obj.kws.get('verbose', False))


if __name__ == '__main__':
    run()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
