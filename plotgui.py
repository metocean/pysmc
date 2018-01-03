#!/usr/bin/env python2.7
# Chris Bunney, Met Office.
# Crown Copyright 2016

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk
from gi.repository import GObject

import matplotlib.pyplot as plt
from matplotlib.backends.backend_gtk3agg import FigureCanvasGTK3Agg as FigureCanvas
import cartopy.crs as ccrs

import netCDF4 as nc
import smc.plotting as smc
import numpy as np

class SMCPlotGui:
    """ PYGtk GUI class wrapper for displaying SMC gridded data """

    def __init__(self, filename=None):
        # initialise some variables:
        self.ncfn = filename
        self.pc = None
        self.cfacs = None
        self.proj = None
        self.src_proj = None

        self.lon1=None
        self.lon2=None
        self.lat1=None
        self.lat2=None

        ###############
        ## GTK setup ##
        ###############

        # create new window
        self.win = Gtk.Window(Gtk.WindowType.TOPLEVEL)

        ######################
        # setup window
        ######################
        self.win.set_border_width(10)
        self.win.set_default_size(600, 400)
        self.win.set_title('SMC Gridded Data Plotter')

        ######################
        # add the GTK canvas:
        ######################
        self.fig = plt.Figure(figsize=(5,4), dpi=100)
        self.canvas = FigureCanvas(self.fig)

        ################
        # Add menu bar #
        ################
        menu_bar = Gtk.MenuBar()
        file_menu = Gtk.Menu()
        open_item = Gtk.MenuItem("Open")
        exit_item = Gtk.MenuItem("Exit")
        file_menu.append(open_item)
        file_menu.append(exit_item)

        open_item.connect("activate", self.load_event)

        root_menu = Gtk.MenuItem("File")
        root_menu.set_submenu(file_menu);
        menu_bar.append(root_menu)


        ###########
        # Controls
        ##########

        # buttons:
        btnPlot = Gtk.Button('Update Plot')
        btnPrev = Gtk.Button('Prev Time')
        btnNext = Gtk.Button('Next Time')

        # Field combo box:
        store = Gtk.ListStore(str,str)
        self.cbox_field = Gtk.ComboBox.new_with_model_and_entry(store)
        cell = Gtk.CellRendererText()
        self.cbox_field.pack_start(cell, True)
        self.cbox_field.add_attribute(cell, 'text', 1)
        self.cbox_field.set_entry_text_column(1)
        store.append(['hs','sig wave heihgt'])


        # Times combo box:
        store = Gtk.ListStore(int,str)
        self.cbox_times = Gtk.ComboBox.new_with_model_and_entry(store)
        cell = Gtk.CellRendererText()
        self.cbox_times.pack_start(cell, True)
        self.cbox_times.add_attribute(cell, 'text', 1)
        self.cbox_times.set_entry_text_column(1)
        #for i in range(1,61):
        #    store.append([i-1, 'T+%03d' % i])

        # Domain combo box:
        store = Gtk.ListStore(str,float,float,float,float)
        self.cbox_domains = Gtk.ComboBox.new_with_model_and_entry(store)
        cell = Gtk.CellRendererText()
        self.cbox_domains.pack_start(cell, True)
        self.cbox_domains.add_attribute(cell, 'text', 0)
        store.append(['Full Domain (could be slow)', -999.9, -999.9, -999.9, -999.9])
        store.append(['UK', 35.0, 70.0, -15.0, 10.0])
        store.append(['South West UK', 49.4, 51.5, -6.7, -1.6])
        store.append(['Mediterranean', 29.5, 46.5, -6.0, 36.5])
        store.append(['North Atlantic', 20.0, 70.0, -90, 30])
        store.append(['West Pacific', -70.0, 70.0, 120, 200])
        store.append(['East Pacific', -70.0, 70.0, -160, -68])
        store.append(['Arabian Gulf', 23.0, 30.5, 47.5, 59.5])
        store.append(['Caspian Sea', 36.0, 47.5, 46.0, 55.5])
        store.append(['Black Sea', 40.5, 47.1, 27.0, 42.0])
        store.append(['Caribbean', 10.0, 27.5, -86.5, -58.5])
        store.append(['South China Sea', -9.5, 24.0, 98.0, 128.0])
        store.append(['Australasia', -48, 0.0, 105.0, 179.0])
        store.append(['New Zealand', -50, -30, 160.0, 182.0])
        self.cbox_domains.set_entry_text_column(0)
        self.cbox_domains.set_active(1)


        # Projections:
        store = Gtk.ListStore(object, str)
        self.cbox_proj = Gtk.ComboBox.new_with_model_and_entry(store)
        cell = Gtk.CellRendererText()
        self.cbox_proj.pack_start(cell, True)
        self.cbox_proj.add_attribute(cell, 'text', 1)
        self.cbox_proj.set_entry_text_column(1)
        store.append([ccrs.PlateCarree(), 'Plate Carree'])
        store.append([ccrs.RotatedPole(pole_latitude=37.5, pole_longitude=177.5),'Euro Rotated Pole'])
        store.append([ccrs.Robinson(), 'Robinson'])
        store.append([ccrs.Mercator(), 'Mercator'])
        store.append([ccrs.Geostationary(), 'Geostationary'])
        store.append([ccrs.PlateCarree(central_longitude=180), 'Plate Carree (central_longitude=180)'])

        self.cbox_proj.set_active(0)

        # coastlines:
        store = Gtk.ListStore(object, str)
        self.cbox_coast = Gtk.ComboBox.new_with_model_and_entry(store)
        cell = Gtk.CellRendererText()
        self.cbox_coast.pack_start(cell, True)
        self.cbox_coast.add_attribute(cell, 'text', 1)
        self.cbox_coast.set_entry_text_column(1)
        store.append([None, 'None'])
        store.append(['10m', 'High res (10m)'])
        store.append(['50m', 'Medium res (50m)'])
        store.append(['110m', 'Low res (110m)'])

        self.cbox_coast.set_active(3)
        self.coast = '110m'

        # lat/lon ranges:
        self.inLat1 = Gtk.Entry();
        self.inLat2 = Gtk.Entry();
        self.inLon1 = Gtk.Entry();
        self.inLon2 = Gtk.Entry();
        self.domain_changed_event(self.cbox_domains) # update with default domain

        # Cell size selection
        cellsbox = Gtk.HBox(homogeneous=False, spacing=5)
        self.chkf1 = Gtk.CheckButton("3km")
        self.chkf2 = Gtk.CheckButton("6km")
        self.chkf3 = Gtk.CheckButton("12km")
        self.chkf4 = Gtk.CheckButton("25km")
        cellsbox.pack_end(self.chkf1, True, True, 0)
        cellsbox.pack_end(self.chkf2, True, True, 0)
        cellsbox.pack_end(self.chkf3, True, True, 0)
        cellsbox.pack_end(self.chkf4, True, True, 0)

        # Colour range box:
        crangebox = Gtk.HBox(homogeneous=False, spacing=5)
        self.cmin = Gtk.Entry()
        self.cmax = Gtk.Entry()
        self.cauto = Gtk.CheckButton('Auto')
        crangebox.pack_start(self.cmin, True, True, 0)
        crangebox.pack_start(self.cmax, True, True, 0)
        crangebox.pack_start(self.cauto, True, True, 0)
        self.cauto.set_active(True)
        self.cmin.set_sensitive(False)
        self.cmax.set_sensitive(False)

        self.cauto.connect('toggled', self.cauto_changed_event)


        ## controls layout
        grid = Gtk.Table(rows=8, columns=3)
        grid.attach(Gtk.Label(label='Field'),     0, 1, 0, 1, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(Gtk.Label(label='Time'),      0, 1, 1, 2, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(Gtk.Label(label='Projection'),0, 1, 2, 3, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(Gtk.Label(label='Coastline'), 0, 1, 3, 4, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(Gtk.Label(label='Domain '),   0, 1, 4, 5, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(Gtk.Label(label='Lat Range'), 0, 1, 5, 6, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(Gtk.Label(label='Lon Range'), 0, 1, 6, 7, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(Gtk.Label(label='Cells'),     0, 1, 7, 8, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(Gtk.Label(label='Colour range'),0, 1, 8, 9, yoptions=Gtk.AttachOptions.SHRINK)

        grid.attach(self.cbox_field,        1, 3, 0, 1, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(self.cbox_times,        1, 3, 1, 2 ,yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(self.cbox_proj,         1, 3, 2, 3 ,yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(self.cbox_coast,        1, 3, 3, 4 ,yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(self.cbox_domains,      1, 3, 4, 5 ,yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(self.inLat1,            1, 2, 5, 6, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(self.inLat2,            2, 3, 5, 6, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(self.inLon1,            1, 2, 6, 7, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(self.inLon2,            2, 3, 6, 7, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(cellsbox,               1, 3, 7, 8, yoptions=Gtk.AttachOptions.SHRINK)
        grid.attach(crangebox,              1, 3, 8, 9, yoptions=Gtk.AttachOptions.SHRINK)


        #grid.attach(btnPlot,                0, 1, 8, 9, yoptions=Gtk.AttachOptions.SHRINK, xoptions=Gtk.AttachOptions.SHRINK)
        # Hbox for plot buttons
        btn_hbox = Gtk.HBox(homogeneous=False, spacing=5)
        btn_hbox.pack_start(btnPrev, True, False, 0)
        btn_hbox.pack_start(btnPlot, True, False, 0)
        btn_hbox.pack_start(btnNext, True, False, 0)


        ## File details text view
        txt = Gtk.TextBuffer()
        txt.set_text('Please load a file')
        self.tv_file_details = Gtk.TextView.new_with_buffer(txt)

        vbox = Gtk.VBox(spacing=15)
        vbox.pack_start(grid, False, True, 0)
        vbox.pack_start(btn_hbox, False, True, 0)
        vbox.pack_end(self.tv_file_details, True, True, 0)

        # plot controls
        from matplotlib.backends.backend_gtk3 import NavigationToolbar2GTK3 as NavigationToolbar
        toolbar = NavigationToolbar(self.canvas, self.win)
        #vbox.pack_end(toolbar, False, False)

        # Top level layout box:
        topbox = Gtk.VBox()
        topbox.pack_start(menu_bar, False, True, 0)

        box = Gtk.HBox(homogeneous=False, spacing=5)
        topbox.pack_end(box, True, True, 0)

        # canvas/toolbar layout
        plotbox = Gtk.VBox(homogeneous=False, spacing=0)
        plotbox.pack_start(self.canvas, True, True, 0)
        plotbox.pack_end(toolbar, False, False, 0)

        box.pack_start(plotbox, True, True, 0)
        box.pack_end(vbox, False, False, 0)
        self.win.add(topbox)


        ###################
        # connect signals:
        ###################
        # destroy/delete:
        self.win.connect("delete_event", self.delete_event)
        self.win.connect("destroy", self.destroy)

        btnPlot.connect("clicked", self.plot_event)
        btnNext.connect("clicked", self.next_time_event)
        btnPrev.connect("clicked", self.prev_time_event)

        self.cbox_domains.connect('changed', self.domain_changed_event)

        # show window
        self.win.show_all()

        #### Load file, if passed in:
        if self.ncfn is not None:
            self.loadfile(self.ncfn)

    def get_cbox_selection(self, combobox, col=0):
        model = combobox.get_model()
        active = combobox.get_active()
        if active < 0:
            return None
        return model[active][col]

    def delete_event(self, widget, event, data=None):
        # returning false after a delete event destroys the widget,
        # returning true means you don't want to destroy the widget.
        return False


    def cauto_changed_event(self, widget, data=None):
        # update the lat/lon boxes:
        active = widget.get_active()
        self.cmin.set_sensitive(not active)
        self.cmax.set_sensitive(not active)

    def domain_changed_event(self, widget, data=None):
        # update the lat/lon boxes:
        model = widget.get_model()
        active = widget.get_active()
        if active < 0:
            return


        for i,o in enumerate([self.inLat1, self.inLat2,
                              self.inLon1, self.inLon2]):
            val = model[active][i+1]
            if val < -900.0:
                val = ''
            else:
                val = str(val)

            o.set_text(val)

    def destroy(self, widget, data=None):
        Gtk.main_quit()

    def load_event(self, widget, data=None):
        # pop up a file selectro dialog and load the file
        if self.selectfile():
            self.loadfile(self.ncfn)

    def selectfile(self):
        dlg = Gtk.FileChooserDialog(title='Select a SMC netCDF file',
                action=Gtk.FileChooserAction.OPEN,
                buttons=(Gtk.STOCK_CANCEL,Gtk.ResponseType.CANCEL,Gtk.STOCK_OPEN,Gtk.ResponseType.OK))
        dlg.set_default_response(Gtk.ResponseType.OK)

        ret = dlg.run()
        if ret == Gtk.ResponseType.OK:
            self.ncfn = dlg.get_filename()

        dlg.destroy()

        return (ret == Gtk.ResponseType.OK)

    def loadfile(self, fn):
        self.d = nc.Dataset(fn, mode='r')

        if self.d.dimensions.has_key('seapoint'):
            seapoint_dim = 'seapoint'
        else:
            seapoint_dim = 'seapoints'

        # populate text bxo with file details
        txt = self.tv_file_details.get_buffer()
        txt.set_text('File details:\n'  +
            '\tNo. sea points: %d\n' % len(self.d.dimensions[seapoint_dim]) +
            '\tNo. times: %d\n' % len(self.d.dimensions['time']))


        # update the fields combobox
        store = self.cbox_field.get_model()
        store.clear()
        store.append(['grid', 'Grid mesh only'])
        for varname,var in self.d.variables.items():
            if var.ndim == 2 and hasattr(var,'long_name'):
                store.append([varname, var.long_name])
            if var.ndim == 3 and hasattr(var,'long_name'):
                store.append(["ens_"+varname, var.long_name + " (ens)"])

        self.cbox_field.set_active(0)


        # find location of uwnd or vwnd and insert a wspd option (this will calculate
        # a wind speed field from uwnd and vwnd):
        i=store.get_iter_first()
        while i is not None:
            varname = store.get_value(i,0)
            if varname in ['uwnd','vwnd']:
                # insert new field here:
                store.insert_before(i, ['wspd_derived', 'wind speed (derived)'])
                break
            i = store.iter_next(i)

        # update times text box:
        store = self.cbox_times.get_model()
        store.clear()
        t = self.d.variables['time']
        t = nc.num2date(t[:], t.units)
        for itime, fctime in enumerate(t):
            store.append([itime, "%s (T+%03d)" % (fctime.strftime("%d/%m/%Y %H:%M"), itime+1)])
        self.cbox_times.set_active(0)

        # get model extents:
        try:
            self.minlat = float(self.d.southernmost_latitude)
            self.maxlat = float(self.d.northernmost_latitude)
            self.minlon = float(self.d.westernmost_longitude)
            self.maxlon = float(self.d.easternmost_longitude)
        except:
            self.minlat = -90
            self.maxlat = 90
            self.minlon = 0
            self.maxlon = 360



        return

    def plot_event(self, widget, data=None):
        self.plotfield()

    def next_time_event(self, widget, data=None):
        # move to next timestep and replot
        itime = self.get_cbox_selection(self.cbox_times)
        if itime is None:
            return

        itime = itime + 1
        ntimes = len(self.cbox_times.get_model())

        if itime >= ntimes:
            itime = 0

        self.cbox_times.set_active(itime);
        self.plot_event(None)

    def prev_time_event(self, widget, data=None):
        # move to next timestep and replot
        itime = self.get_cbox_selection(self.cbox_times)
        if itime is None:
            return

        itime = itime - 1
        ntimes = len(self.cbox_times.get_model())

        if itime < 0:
            itime = ntimes - 1

        self.cbox_times.set_active(itime);
        self.plot_event(None)


    def getfloat(self, string):
        try:
            return float(string)
        except ValueError:
            return None

    def plotfield(self):

        # get reference to details text view and clear it
        txt = self.tv_file_details.get_buffer()
        start_iter = txt.get_start_iter()
        end_iter = txt.get_end_iter()
        txt.delete(start_iter, end_iter)
        end_iter = txt.get_end_iter()


        # get required field, time, etc:
        var = self.get_cbox_selection(self.cbox_field)
        if var.startswith('ens_'):
            var = var[4:]
            ens = True
        else:
            ens = False

        itime = self.get_cbox_selection(self.cbox_times)


        # get lat/lon ranges
        lat1 = self.getfloat(self.inLat1.get_text())
        lat2 = self.getfloat(self.inLat2.get_text())
        lon1 = self.getfloat(self.inLon1.get_text())
        lon2 = self.getfloat(self.inLon2.get_text())

        # get celfacs:
        cfacs=[]
        if self.chkf1.get_active():
            cfacs.append(1)
        if self.chkf2.get_active():
            cfacs.append(2)
        if self.chkf3.get_active():
            cfacs.append(4)
        if self.chkf4.get_active():
            cfacs.append(8)

        if len(cfacs) == 0:
            cfacs = None

        if cfacs != self.cfacs:
            self.pc = None
            self.cfacs = cfacs


        if(lat1 != self.lat1 or lat2 != self.lat2 or
                lon1 != self.lon1 or lon2 != self.lon2 ):

            self.pc = None
            self.lat1 = lat1
            self.lat2 = lat2
            self.lon1 = lon1
            self.lon2 = lon2

        proj = self.get_cbox_selection(self.cbox_proj)
        if self.proj != proj:
            self.proj = proj
            self.pc = None
            self.fig.clf()
            self.ax = self.fig.add_subplot(111, projection=self.proj)
            sm = plt.cm.ScalarMappable(norm=plt.Normalize())
            sm._A=[]
            self.cbar = plt.colorbar(sm, ax=self.ax, orientation='horizontal', fraction=0.05, shrink=0.8, pad=0.04)

        coast = self.get_cbox_selection(self.cbox_coast)
        if self.coast != coast:
            newcoast = True
            self.coast = coast
        else:
            newcoast = False

        ax = self.ax

        # determin extents:
        if isinstance(proj, ccrs.Geostationary) or (
                lon1 is None and lon2 is None and lat1 is None and lat2 is None):
            ax.set_global()
            global_extent = True
        else:
            if lon1 is None:
                lon1 = self.minlon
            if lon2 is None:
                lon2 = self.maxlon
            if lat1 is None:
                lat1 = self.minlat
            if lat2 is None:
                lat2 = self.maxlat
            global_extent = False


        ## check source projection of requested variable:
        chkvar = var
        if var == 'wspd_derived': chkvar = 'uwnd'
        if var == 'grid': chkvar = 'hs'
        if hasattr(self.d.variables[chkvar], "grid_mapping"):
            # get rotated pole projection:
            mapname = self.d.variables[chkvar].grid_mapping
            gmap = self.d.variables[mapname]
            plat = gmap.grid_north_pole_latitude
            plon = gmap.grid_north_pole_longitude
            src_proj = ccrs.RotatedPole(pole_latitude=plat, pole_longitude=plon)
        else:
            src_proj = None

        if self.src_proj != src_proj:
            # source projection changed need to update patches:
            self.pc = None
            self.src_proj = src_proj


        # If patch collection not yet calculated, then build it now and add it to
        # the axis along with the coastlines:
        if self.pc is None:
            txt.insert(end_iter, "Domain has changed - recalculating patches...")
            while Gtk.events_pending(): Gtk.main_iteration_do(True)

            self.pc = smc.generate_patch_collection(self.d,
                    lon1=self.lon1, lon2=self.lon2,
                    lat1=self.lat1, lat2=self.lat2,
                    target_crs=self.proj, cfacs=cfacs,
                    source_crs=self.src_proj)

            end_iter = txt.get_end_iter()
            txt.insert(end_iter, "Done\n")
            while Gtk.events_pending(): Gtk.main_iteration_do(True)

            end_iter = txt.get_end_iter()
            txt.insert(end_iter, "Adding patches to axes\n");
            print("Adding patches to axes\n")
            while Gtk.events_pending(): Gtk.main_iteration_do(True)

            ax.cla()
            self.cm = ax.add_collection(self.pc.pcol)

            if not global_extent:
                ax.set_extent((lon1, lon2, lat1, lat2), crs=ccrs.PlateCarree())

            if self.coast is not None:
                end_iter = txt.get_end_iter()
                txt.insert(end_iter, "Adding coastline\n");
                while Gtk.events_pending(): Gtk.main_iteration_do(True)
                ax.coastlines(resolution=self.coast)

        elif newcoast:
            # coastline changed; clear axis and redraw:
            ax.cla()

            end_iter = txt.get_end_iter()
            txt.insert(end_iter, "Adding patches to axes\n");
            while Gtk.events_pending(): Gtk.main_iteration_do(True)
            self.cm = ax.add_collection(self.pc.pcol)

            if not global_extent:
                ax.set_extent((lon1, lon2, lat1, lat2), crs=ccrs.PlateCarree())

            if self.coast is not None:
                end_iter = txt.get_end_iter()
                txt.insert(end_iter, "Adding coastline\n");
                while Gtk.events_pending(): Gtk.main_iteration_do(True)
                ax.coastlines(resolution=self.coast)


        # update the patch face colours with the new data:
        end_iter = txt.get_end_iter()
        txt.insert(end_iter, "Updating patches...\n");
        while Gtk.events_pending(): Gtk.main_iteration_do(True)

        if var == 'grid':
            # special case - just plot grid mesh:
            self.pc.pcol.set_facecolor('#EAEAEA')
            self.pc.pcol.set_edgecolor('#5E5E5E')
            self.pc.pcol.set_linewidth(0.5)
            cmin = 0
            cmax = 1
        else:
            # load field and set face colours
            if var == 'wspd_derived':
                end_iter = txt.get_end_iter()
                txt.insert(end_iter, "Deriving wind speed from uwnd and vwnd fields\n")
                fld = self.d.variables['uwnd']
                fld2 = self.d.variables['vwnd']
                dat = np.sqrt(
                        np.power(fld[itime,:][self.pc.mask], 2) +
                        np.power(fld2[itime,:][self.pc.mask], 2) )
            else:
                fld = self.d.variables[var]
                if ens:
                    dat = fld[0,itime,:][self.pc.mask]  # just get member 0
                else:
                    dat = fld[itime,:][self.pc.mask]

            if self.cauto.get_active():
                cmin = dat.min()
                cmax = dat.max()
                self.cmin.set_text("%.2f" % cmin)
                self.cmax.set_text("%.2f" % cmax)
            else:
                cmin = self.getfloat(self.cmin.get_text())
                cmax = self.getfloat(self.cmax.get_text())

            clrs = smc.generate_color_array(dat, cmin=cmin, cmax=cmax)
            self.pc.pcol.set_facecolor(clrs)
            self.pc.pcol.set_edgecolor(clrs)
            self.pc.pcol.set_linewidth(0.5)


        # update colorbar mappable with new range:
        self.cbar.mappable.set_clim(cmin, cmax)

        ## update title
        fldname = self.get_cbox_selection(self.cbox_field, 1)
        time = self.get_cbox_selection(self.cbox_times, 1)
        ax.set_title("%s @ %s" % (fldname, time), fontsize=12)

        # update canvas
        self.canvas.draw()
        end_iter = txt.get_end_iter()
        txt.insert(end_iter, "Done.\n");
        while Gtk.events_pending(): Gtk.main_iteration_do(True)

    def main(self):
        Gtk.main()


if __name__ == '__main__':
    #gui = SMCPlotGui('/home/h06/frey/wvgu_out_grd.nc1')
    import sys
    fn = None
    if len(sys.argv) > 1:
        fn = sys.argv[1]
    gui = SMCPlotGui(fn)
    gui.main()
