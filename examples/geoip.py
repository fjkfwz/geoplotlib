# coding: utf-8

import csv
import random
import numpy as np
import geoplotlib

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import NullFormatter
from matplotlib.mlab import csv2rec
from matplotlib.cbook import get_sample_data
from geoplotlib.utils import BoundingBox
from mapapi import baidu
from collections import defaultdict
from geoplotlib.utils import DataAccessObject


class GetGip(object):
    def __init__(self):
        self.csvfile = 'Y公佈20151213.csv'
        self.gipfile = 'data1.csv'
        self.map_api = baidu.MapApi(['859b06c7ebdaa6301dafbeea74a57677'])

    def get_and_save(self):
        wf = open(self.gipfile, "wb")
        with open(self.csvfile) as f:
            reader = csv.reader(f, delimiter=',', quotechar='"')
            writer = csv.writer(wf)
            for row in reader:
                if not row:
                    break
                if row[0] == "近祖姓":
                    continue
                p = row[2]
                l = row[3]

                location = self.map_api.location_api.get_location_by_address(unicode(l, "utf-8"), unicode(p, "utf-8"))
                if not location:
                    continue
                # print type(location)
                # dict = simplejson.loads(location)
                # print location
                # print location['lat'],location['lng']
                writer.writerow(row + [location['lat'], location['lng']])


class CaculateDateAndSave(object):
    def __init__(self, source_file, output_filename):
        self.source_file = source_file
        self.normalization_lon_haplogroup_dict = defaultdict(lambda: defaultdict(int))
        self.normalization_lat_haplogroup_dict = defaultdict(lambda: defaultdict(int))
        self.lat_haplogroup_dict = defaultdict(lambda: defaultdict(int))
        self.lon_haplogroup_dict = defaultdict(lambda: defaultdict(int))
        self.ouput_filename = output_filename

    def load_data(self):
        with open(self.source_file) as f:
            reader = csv.reader(f, delimiter=',', quotechar='"')
            for row in reader:
                yield row

    def get_haplogroup_amount(self):
        for row in self.load_data():
            # lon = ('lon', row[-1])
            # lat = ('lat', row[-2])
            lon = row[-1]
            lat = row[-2]
            loca = row[3]
            haplogroup_classfy = str(row[4] or row[6])
            if not haplogroup_classfy:
                continue
            d_haplogroup_classfy = [haplogroup_classfy[0:i + 1] for i in range(len(haplogroup_classfy))]

            for haplo in d_haplogroup_classfy:
                self.lat_haplogroup_dict[haplo][lat] += 1
                self.lon_haplogroup_dict[haplo][lon] += 1

    def write_to_file(self):
        self.get_haplogroup_amount()
        self.normalization_by_lat()
        self.normalization_by_lon()
        lat_f = open("lat_%s" % self.ouput_filename, "wb")
        lat_writer = csv.writer(lat_f)
        for (haplo, lat_list) in self.normalization_lat_haplogroup_dict.items():
            for lat, amount in lat_list.items():
                lat_writer.writerow([haplo, lat, amount])
        lat_f.close()

        lon_f = open("lon_%s" % self.ouput_filename, "wb")
        lon_writer = csv.writer(lon_f)
        for (haplo, lon_list) in self.normalization_lon_haplogroup_dict.items():
            for lon, amount in lon_list.items():
                lon_writer.writerow([haplo, lon, amount])
        lon_f.close()

    def normalization_by_lat(self):
        cnt_dict = defaultdict(int)
        for (haplo, lat_list) in self.lat_haplogroup_dict.items():
            for (lat, number) in lat_list.items():
                if len(haplo) > 1:
                    continue
                cnt_dict[lat] += number

        max_cnt = max(cnt_dict.values())
        print max_cnt

        time = dict((lat, float(cnt) / max_cnt) for lat, cnt in cnt_dict.items())

        # print time
        for (haplo, lat_list) in self.lat_haplogroup_dict.items():
            for (lat, number) in lat_list.items():
                self.normalization_lat_haplogroup_dict[haplo][lat] = float(number) / time[lat]

    def normalization_by_lon(self):
        cnt_dict = defaultdict(int)
        for (haplo, lon_list) in self.lon_haplogroup_dict.items():
            for (lon, number) in lon_list.items():
                if len(haplo) > 1:
                    continue
                cnt_dict[lon] += number

        max_cnt = max(cnt_dict.values())
        print max_cnt

        time = dict((lon, float(cnt) / max_cnt) for lon, cnt in cnt_dict.items())

        # print time
        for (haplo, lon_list) in self.lon_haplogroup_dict.items():
            for (lon, number) in lon_list.items():
                self.normalization_lon_haplogroup_dict[haplo][lon] = float(number) / time[lon]

    def get_dealed_lat_date(self):
        lat_list = self.normalization_by_lat().keys()
        lon_list = self.normalization_by_lon().keys()
        min_lat = min(lat_list)
        min_lon = min(lon_list)
        max_lon = max(lon_list)
        max_lat = max(lat_list)
        lat_interval = max_lat - min_lat
        lon_interval = max_lon - min_lon

    def return_all_lat(self):
        cnt_dict = defaultdict(lambda: defaultdict(int))
        for (haplo, haplo_list) in self.lat_haplogroup_dict.items():
            for (lat, distributed) in haplo_list.items():
                cnt_dict[lat][haplo] = distributed
        return cnt_dict

    def return_all_lon(self):
        cnt_dict = defaultdict(lambda: defaultdict(int))
        for (haplo, haplo_list) in self.lon_haplogroup_dict.items():
            for (lon, distributed) in haplo_list.items():
                cnt_dict[lon][haplo] = distributed
        return cnt_dict

    def stacked_draw(self):
        N = 5
        menMeans = (20, 35, 30, 35, 27)
        womenMeans = (25, 32, 34, 20, 25)
        menStd = (2, 3, 4, 1, 2)
        womenStd = (3, 5, 2, 3, 3)
        ind = np.arange(N)    # the x locations for the groups
        width = 0.35       # the width of the bars: can also be len(x) sequence

        p1 = plt.bar(ind, menMeans, width, color='r', yerr=menStd)
        p2 = plt.bar(ind, womenMeans, width, color='y',
                     bottom=menMeans, yerr=womenStd)

        plt.ylabel('Scores')
        plt.title('Scores by group and gender')
        plt.xticks(ind + width/2., ('G1', 'G2', 'G3', 'G4', 'G5'))
        plt.yticks(np.arange(0, 81, 10))
        plt.legend((p1[0], p2[0]), ('Men', 'Women'))

        plt.show()

    def draw(self):
        self.get_haplogroup_amount()
        self.normalization_by_lat()
        self.normalization_by_lon()
        # the random data
        cmap_list = ['Blues', 'BuGn', 'BuPu',
                     'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
                     'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
                     'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']

        colors = ['c', 'y', 'm', 'r', 'mediumseagreen', 'crimson', 'orange', 'teal', 'mediumpurple', 'indigo']

        g_x = map(float, self.return_all_lon().keys())
        medium_x = min(g_x) + (max(g_x) - min(g_x)) / 2
        g_x = [p - medium_x for p in g_x]
        g_y = map(float, self.return_all_lat().keys())
        medium_y = min(g_y) + (max(g_y) - min(g_y)) / 2
        g_y = [p - medium_y for p in g_y]

        nullfmt = NullFormatter()  # no labels

        # definitions for the axes
        left, width = 0.05, 0.65
        bottom, height = 0.05, 0.65
        bottom_h = left_h = left + width + 0.02

        rect_scatter = [left, bottom, width, height]
        rect_histx = [left, bottom_h, width, 0.22]
        rect_histy = [left_h, bottom, 0.22, height]
        rect_label = [left_h, bottom_h, 0.22, 0.22]

        # start with a rectangular Figure
        plt.figure(1, figsize=(12, 10))

        axScatter = plt.axes(rect_scatter)
        axHistx = plt.axes(rect_histx)
        axHisty = plt.axes(rect_histy)
        axLable = plt.axes(rect_label, frameon=False)

        # no labels
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)
        axLable.yaxis.set_major_formatter(nullfmt)
        axLable.xaxis.set_major_formatter(nullfmt)

        # now determine nice limits by hand:
        binwidth = 1
        xmax = np.max(np.fabs(g_x))
        xlim = ((int(xmax / binwidth) + 1) * binwidth) * 1.1

        ymax = np.max(np.fabs(g_y))
        ylim = ((int(ymax / binwidth) + 1) * binwidth) * 1.1

        axScatter.set_xlim((-xlim, xlim))
        axScatter.set_ylim((-ylim, ylim))

        xbins = np.arange(-xlim, xlim + binwidth, binwidth)
        ybins = np.arange(-ylim, ylim + binwidth, binwidth)
        patch_list = list()
        for haplo, lon_dict in self.normalization_lon_haplogroup_dict.items():
            if len(haplo) > 1 or haplo == "?":
                continue
            color = colors.pop()
            patch_list.append(mpatches.Patch(color=color, label=haplo, alpha=0.5))
            lat_dict = self.normalization_lat_haplogroup_dict[haplo]
            x = map(float, lon_dict.keys())
            x = [p - medium_x for p in x]
            y = map(float, lat_dict.keys())
            y = [p - medium_y for p in y]
            s = lon_dict.values()
            # the scatter plot:
            for p in xrange(len(x)):
                # print x[p], y[p], s[p]
                axScatter.scatter(x[p], y[p], alpha=0.5, s=s[p] * 20, cmap='Greys', color=color)

            x = list()
            for lon, number in self.normalization_lon_haplogroup_dict[haplo].items():
                x += [float(lon) - medium_x] * int(number)

            axHistx.hist(x, bins=xbins, color=color, alpha=0.5)

            y = list()
            for lat, number in self.normalization_lat_haplogroup_dict[haplo].items():
                y += [float(lat) - medium_y] * int(number)

            axHisty.hist(y, bins=ybins, color=color, orientation='horizontal', alpha=0.5)

        axHistx.set_xlim(axScatter.get_xlim(), auto=True)
        axHisty.set_ylim(axScatter.get_ylim(), auto=True)
        axLable.legend(handles=patch_list, mode="expand", ncol=2, shadow=True)
        axHistx.set_title('Distribution Of Haplogroup %s' % "", fontsize=21, ha='center', va='bottom', style="oblique")
        plt.draw()
        plt.savefig('%s_distributed.png' % "",  bbox_inches='tight')
        plt.show()
        plt.close()


class DrawMaps(object):
    def __init__(self, data, cmap, r, output_filename):
        self.data = data
        self.cmap = cmap
        self.r = r
        self.ouput_filename = output_filename

    def get_dao_object(self):
        values = defaultdict(list)
        for (k, v) in self.data:
            values[k].append(v)
        npvalues = {k: np.array(values[k]) for k in values.keys()}
        for k in npvalues.keys():
            for datatype in [np.int, np.float]:
                try:
                    npvalues[k][:1].astype(datatype)
                    npvalues[k] = npvalues[k].astype(datatype)
                    break
                except:
                    pass
        dao = DataAccessObject(npvalues)
        return dao

    def drawmap(self):
        """
        Multiple examples of kernel density estimation visualization
        """
        data = self.get_dao_object()
        # geoplotlib.kde(data, bw=5, cut_below=1e-4)

        # lowering clip_above changes the max value in the color scale
        # geoplotlib.kde(data, bw=5, cut_below=1e-4, clip_above=.1)

        # different bandwidths
        geoplotlib.kde(data, bw=20, cmap=self.cmap, cut_below=1e-4)
        # geoplotlib.kde(data, bw=2, cmap='PuBuGn', cut_below=1e-4)

        # linear colorscale
        # geoplotlib.kde(data, bw=5, cmap='jet', cut_below=1e-4, scaling='lin')

        geoplotlib.set_bbox(BoundingBox.from_nominatim('CHINA'))

        geoplotlib.savefig(self.ouput_filename)


class SpiltData(object):
    def __init__(self, source_file):
        self.source_file = source_file
        self.name_dict = defaultdict(list)
        self.haplogroup_dict = defaultdict(list)
        self.snp_dict = defaultdict(list)
        self.province_dict = defaultdict(lambda: defaultdict(list))
        self.location_dict = defaultdict(lambda: defaultdict(list))
        self.p_haplogroup_dict = defaultdict(list)
        self.d_province_dict = defaultdict(lambda: defaultdict(list))
        self.normalization_location_dict = defaultdict(lambda: defaultdict(list))
        self.normalization_province_dict = defaultdict(lambda: defaultdict(list))
        self.normalization_d_province_dict = defaultdict(lambda: defaultdict(list))

    def load_data(self):
        with open(self.source_file) as f:
            reader = csv.reader(f, delimiter=',', quotechar='"')
            for row in reader:
                yield row

    def spilt_data(self):
        for row in self.load_data():
            lon = ('lon', row[-1])
            lat = ('lat', row[-2])
            loca = row[3]
            haplogroup_classfy = str(row[4] or row[6])
            if not haplogroup_classfy:
                continue
            p_haplogroup_classfy = haplogroup_classfy[0]

            d_haplogroup_classfy = [haplogroup_classfy[0:i + 1] for i in range(len(haplogroup_classfy))]

            if not haplogroup_classfy:
                continue
            for haplo in d_haplogroup_classfy:
                self.d_province_dict[loca][haplo].append(lat)
                self.d_province_dict[loca][haplo].append(lon)

            self.location_dict[loca][p_haplogroup_classfy].append(lat)
            self.location_dict[loca][p_haplogroup_classfy].append(lon)

            # self.province_dict[loca][p_haplogroup_classfy].append(lat)
            # self.province_dict[loca][p_haplogroup_classfy].append(lon)

            self.province_dict[loca][haplogroup_classfy].append(lat)
            self.province_dict[loca][haplogroup_classfy].append(lon)

            self.p_haplogroup_dict[p_haplogroup_classfy].append(lat)
            self.p_haplogroup_dict[p_haplogroup_classfy].append(lon)

            self.haplogroup_dict[haplogroup_classfy].append(lat)
            self.haplogroup_dict[haplogroup_classfy].append(lon)

            # 不同民族的单倍型分布

    def normalization_by_location(self):
        cnt_dict = defaultdict(int)
        for (loca, haplogroup) in self.location_dict.items():
            for (haplo, distributed) in haplogroup.items():
                cnt_dict[loca] += len(distributed)

        max_cnt = max(cnt_dict.values())

        time = dict((loca, float(cnt) / max_cnt) for loca, cnt in cnt_dict.items())

        # print time
        for (loca, haplogroup) in self.location_dict.items():
            for (haplo, distributed) in haplogroup.items():
                self.normalization_location_dict[loca][haplo] = distributed * int(float(1) / time[loca])

    def normalization_by_province(self):
        cnt_dict = defaultdict(int)
        for (loca, haplogroup) in self.province_dict.items():
            for (haplo, distributed) in haplogroup.items():
                cnt_dict[loca] += len(distributed)

        max_cnt = max(cnt_dict.values())
        print max_cnt

        time = dict((loca, float(cnt) / max_cnt) for loca, cnt in cnt_dict.items())

        # print time
        for (loca, haplogroup) in self.province_dict.items():
            for (haplo, distributed) in haplogroup.items():
                self.normalization_province_dict[loca][haplo] = distributed * int(float(1) / time[loca])

    def normalization_by_d_province(self):
        cnt_dict = defaultdict(int)
        for (loca, haplogroup) in self.d_province_dict.items():
            for (haplo, distributed) in haplogroup.items():
                cnt_dict[loca] += len(distributed)

        max_cnt = max(cnt_dict.values())
        print max_cnt

        time = dict((loca, float(cnt) / max_cnt) for loca, cnt in cnt_dict.items())

        # print time
        for (loca, haplogroup) in self.d_province_dict.items():
            for (haplo, distributed) in haplogroup.items():
                self.normalization_d_province_dict[loca][haplo] = distributed * int(float(1) / time[loca])

    def draw(self):
        cmap_list = ['Blues', 'BuGn', 'BuPu',
                     'GnBu', 'Greens', 'Greys', 'Oranges', 'OrRd',
                     'PuBu', 'PuBuGn', 'PuRd', 'Purples', 'RdPu',
                     'Reds', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd']

        # 绘制单倍群分布
        # for filename, item in self.p_haplogroup_dict.items():
        #     DrawMaps(item, cmap_list[1], None, filename).drawmap()

        # 绘制normalization后的单倍群分布
        p_haplogroup_dict = defaultdict(list)
        # for (loca, haplogroup) in self.normalization_location_dict.items():
        #     for (haplo, distributed) in haplogroup.items():
        #         p_haplogroup_dict[haplo] += distributed
        # for (filename, item) in p_haplogroup_dict.items():
        #     DrawMaps(item, cmap_list[1], None, filename).drawmap()

        for (loca, haplogroup) in self.normalization_d_province_dict.items():
            for (haplo, distributed) in haplogroup.items():
                p_haplogroup_dict[haplo] += distributed
        for (filename, item) in p_haplogroup_dict.items():
            DrawMaps(item, random.choice(cmap_list), None, filename).drawmap()

    def inter(self):
        self.spilt_data()
        self.normalization_by_d_province()
        # self.draw()

    def draw_lat_haplogroup_distributed(self):
        fname = get_sample_data('percent_bachelors_degrees_women_usa.csv')
        gender_degree_data = csv2rec(fname)

        # These are the colors that will be used in the plot
        color_sequence = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c',
                          '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5',
                          '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f',
                          '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']

        # You typically want your plot to be ~1.33x wider than tall. This plot
        # is a rare exception because of the number of lines being plotted on it.
        # Common sizes: (10, 7.5) and (12, 9)
        fig, ax = plt.subplots(1, 1, figsize=(12, 14))

        # Remove the plot frame lines. They are unnecessary here.
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        # Ensure that the axis ticks only show up on the bottom and left of the plot.
        # Ticks on the right and top of the plot are generally unnecessary.
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()

        # Limit the range of the plot to only where the data is.
        # Avoid unnecessary whitespace.
        plt.xlim(1968.5, 2011.1)
        plt.ylim(-0.25, 90)

        # Make sure your axis ticks are large enough to be easily read.
        # You don't want your viewers squinting to read your plot.
        plt.xticks(range(1970, 2011, 10), fontsize=14)
        plt.yticks(range(0, 91, 10), ['{0}%'.format(x)
                                      for x in range(0, 91, 10)], fontsize=14)

        # Provide tick lines across the plot to help your viewers trace along
        # the axis ticks. Make sure that the lines are light and small so they
        # don't obscure the primary data lines.
        for y in range(10, 91, 10):
            plt.plot(range(1969, 2012), [y] * len(range(1969, 2012)), '--',
                     lw=0.5, color='black', alpha=0.3)

        # Remove the tick marks; they are unnecessary with the tick lines we just
        # plotted.
        plt.tick_params(axis='both', which='both', bottom='off', top='off',
                        labelbottom='on', left='off', right='off', labelleft='on')

        # Now that the plot is prepared, it's time to actually plot the data!
        # Note that I plotted the majors in order of the highest % in the final year.
        majors = ['Health Professions', 'Public Administration', 'Education',
                  'Psychology', 'Foreign Languages', 'English',
                  'Communications\nand Journalism', 'Art and Performance', 'Biology',
                  'Agriculture', 'Social Sciences and History', 'Business',
                  'Math and Statistics', 'Architecture', 'Physical Sciences',
                  'Computer Science', 'Engineering']

        y_offsets = {'Foreign Languages': 0.5, 'English': -0.5,
                     'Communications\nand Journalism': 0.75,
                     'Art and Performance': -0.25, 'Agriculture': 1.25,
                     'Social Sciences and History': 0.25, 'Business': -0.75,
                     'Math and Statistics': 0.75, 'Architecture': -0.75,
                     'Computer Science': 0.75, 'Engineering': -0.25}

        for rank, column in enumerate(majors):
            # Plot each line separately with its own color.
            column_rec_name = column.replace('\n', '_').replace(' ', '_').lower()

            line = plt.plot(gender_degree_data.year,
                            gender_degree_data[column_rec_name],
                            lw=2.5,
                            color=color_sequence[rank])

            # Add a text label to the right end of every line. Most of the code below
            # is adding specific offsets y position because some labels overlapped.
            y_pos = gender_degree_data[column_rec_name][-1] - 0.5

            if column in y_offsets:
                y_pos += y_offsets[column]

            # Again, make sure that all labels are large enough to be easily read
            # by the viewer.
            plt.text(2011.5, y_pos, column, fontsize=14, color=color_sequence[rank])

        # Make the title big enough so it spans the entire plot, but don't make it
        # so big that it requires two lines to show.

        # Note that if the title is descriptive enough, it is unnecessary to include
        # axis labels; they are self-evident, in this plot's case.
        plt.title('Percentage of Bachelor\'s degrees conferred to women in '
                  'the U.S.A. by major (1970-2011)\n', fontsize=18, ha='center')

        # Finally, save the figure as a PNG.
        # You can also save it as a PDF, JPEG, etc.
        # Just change the file extension in this call.
        plt.savefig('percent-bachelors-degrees-women-usa.png', bbox_inches='tight')

    def draw_lon_haplogroup_distributed(self):
        pass


if __name__ == "__main__":
    # GetGip().get_and_save()
    # SpiltData('data1.csv').inter()
    CaculateDateAndSave('data1.csv', 'haplo_quality').draw()
