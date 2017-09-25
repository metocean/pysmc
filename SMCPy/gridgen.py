import subprocess as subp
import sys
import os
import logging as log
import docker
from shutil import copy

log.basicConfig(level=log.INFO)

def create_grid(id,latmin,latmax,lonmin,lonmax,dlat,dlon,bathy='etopo1',
        coastalres='full',isglobal=0,lake_tol=None,
        outdir='./out', overwrite=False, testmode=0,
        userpoly_file=None, shift_userpolys=False,
        matcmd='create_grid_smcbase', mdocker=True):
    """
    latmin          :; Minimum latatide
    lonmin          :; Minimum longitude
    latmax          :; Maximum latatide
    lonmax          :; Maximum longitude
    dlat            :: latitude increment
    dlon            :: longitude increment
    bathy           :: bathy source
    coastalres      :: coastal resolution for cartopy coastal polygons
    isglobal        :: grid is global (0 - False, 1 - True)
    lake_tol        :: lake_tol (see gridgen docs)
    outdir          :: output directory
    overwrite       :: overwrite existing files
    testmode        :: don't consider polyons in mas calculation (speed at the
                       cost of accuracy)
    userpoly_file   :: User polygon file
    shift_userpolys :: Shift userpolys to 0,360 range
    mdocker         :: Run matlab command in matlab docker
    """

    metafile=os.path.join(outdir,id+'.meta')
    if os.path.isfile(metafile):
        log.info("File %s exists" % metafile)
        if overwrite:
            log.info("Overwriting files")
        else:
            return
    else:
        if not os.path.isdir(outdir):
            os.makedirs(outdir)

    # If not specified, set default lake_tol to -1 (largest body of water only)
    if not lake_tol:
        lake_tol = -1
    if not userpoly_file:
        fname_poly = 0
    else:
        fname_poly = userpoly_file

    # Gridgen likes lons in range 0, 360 (this may need to be looked)
    # lonmin%=360.
    # lonmax%=360.

    cmd = ['matlab','-nodesktop', '-nosplash', '-r']
    matargs = []
    matargs.append(matcmd)
    matargs.append(id)
    matargs.append(str(latmin))
    matargs.append(str(latmax))
    matargs.append(str(lonmin))
    matargs.append(str(lonmax))
    matargs.append(str(dlat))
    matargs.append(str(dlon))
    matargs.append(bathy)
    matargs.append(coastalres)
    matargs.append(str(isglobal))
    matargs.append(str(lake_tol))
    if mdocker:
        matargs.append('/tmp')
    else:
        matargs.append(outdir)
    matargs.append(str(testmode))
    matargs.append(str(fname_poly))
    matargs.append(str(int(shift_userpolys)))
    cmd.append(' '.join(matargs))
    log.info("Running " + ' '.join(cmd))
    if mdocker:
        DockerRun(cmd)
        copy('/tmp/'+id+'.mat', outdir)
    else:
        subp.call(cmd)

class DockerRun(object):

    def __init__(self,
                 cmd,
                 headip='localhost',
                 dopts={'image': 'metocean/matlab',
                        'user': 'metocean',
                        'volumes': ['/flush:/flush',
                                    '/source:/source',
                                    '/source/pysmc/SMCPy/matlab:/home/metocean/Documents/MATLAB',
                                    '/tmp:/tmp',
                                    ],
                        'mac_address': '02:42:ac:11:00:05',
                        'working_dir': '/tmp',
                        },
                 logger=log,
                 **kwargs):
        self.headclient = docker.DockerClient(base_url='tcp://%s:2375' % headip)
        self.dopts = dopts
        self.cmd = cmd
        self.logger=logger
        self.run()

    def run_cmd(self):
        self.logger.info("Running %s" % self.cmd)
        output = self.head.exec_run(self.cmd, stream=True, user='metocean')
        for val in output:
            self.logger.info(val.strip())

    def pull(self):
        self.logger.info("Pulling containers...")
        self.logger.info("  ...head")
        self.head = self.headclient.images.pull(self.dopts['image'])

    def down(self):
        """Need to kill the ssh daemons as these are sudo processes that hang
        the docker stop"""
        self.logger.info("Bringing down containers...")
        self.head.stop()
        self.head.remove(force=True)

    def up(self):
        self.logger.info("Bringing up containers...")
        self.logger.info("  ...head")
        self.head = self.headclient.containers.run(detach=True, command='sleep infinity', **self.dopts)

    def run(self):
        try:
            self.pull()
            self.up()
            self.run_cmd()
        except Exception as e:
            raise e
        finally:
            self.down()
