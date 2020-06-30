from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph
import logging

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Alfa', anchor='alfa',
                                            href="https://github.com/biocompibens/ALFA",
                                            info="is an example analysis module used for writing documentation.")

        self.alfa_data = dict()
        self.alfa_data['tempp'] = {}

        self.categories = dict()
        self.categories['tempp'] = {}

        self.biotypes = dict()
        self.biotypes['tempp'] = {}

        categories_temp = dict()
        categories_temp['tempp'] = {}

        biotypes_temp = dict()
        biotypes_temp['tempp'] = {}

        # Find and load any ALFA reports
        number = 0

        filesDone = []

        for f in self.find_log_files('alfa'):

            alfaFilename = self.getFilename(f)

            if alfaFilename not in filesDone:
                number = number + 1

                filesDone = filesDone + [alfaFilename]

                categories_temp, biotypes_temp = self.parse_alfa_logs(
                    f, alfaFilename)
                self.updateDictValues(
                    categories_temp, biotypes_temp, alfaFilename)
            #print(self.categories, self.biotypes)

        if number == 0:
            raise UserWarning
        # self.calculatePercentage()
        #print(self.categories, self.biotypes)
        self.print_alfa_charts()

    def getFilename(self, f):
        fullFilename = f['fn']
        part2 = fullFilename.find(".ALFA_feature_counts.tsv")
        part1 = fullFilename[0:part2]
        return part1

    def parse_alfa_logs(self, f, alfaFilename):
        j = 0
        word = []
        words = []
        listToStr1 = []
        listToStr = []
        for l in f['f']:
            if(l != '\t' and l != '\n'):
                # print(l)
                word.append(l)
            else:
                words.append(word)
                word = []
        for k in words:
            listToStr1.append(''.join([str(elem) for elem in k]))
        for k in listToStr1:
            listToStr.append(''.join([str(elem) for elem in k]))
        # print(listToStr)

        catdict = dict()
        catdict['tempp'] = {}

        biodict = dict()
        biodict['tempp'] = {}

        for i_num in range(len(listToStr)):
            if i_num < 3:
                continue
            elif i_num % 3 == 0:
                cat, bio = listToStr[i_num].split(',')

                catdict.setdefault(alfaFilename, {})
                biodict.setdefault(alfaFilename, {})

                catdict[alfaFilename].setdefault(cat, 0)
                biodict[alfaFilename].setdefault(bio, 0)

                catdict[alfaFilename][cat] = catdict[alfaFilename
                                                     ][cat] + float(listToStr[i_num + 1])
                biodict[alfaFilename][bio] = biodict[alfaFilename
                                                     ][bio] + float(listToStr[i_num + 1])

        return catdict, biodict

    def updateDictValues(self, categories_temp, biotypes_temp, alfaFilename):
        if 'tempp' in self.categories:
            del self.categories['tempp']
        if 'tempp' in self.biotypes:
            del self.biotypes['tempp']
        del categories_temp['tempp']
        del biotypes_temp['tempp']

        # Combining both the dictionaries
        for i in categories_temp.keys():
            found = 0
            for j in self.categories.keys():
                if(i == j):
                    found = 1
                    self.categories[j].update(categories_temp[i])
            if(found == 0):
                self.categories[i] = categories_temp[i]

        for i in biotypes_temp.keys():
            found = 0
            for j in self.biotypes.keys():
                if(i == j):
                    found = 1
                    self.biotypes[j].update(biotypes_temp[i])
            if(found == 0):
                self.biotypes[i] = biotypes_temp[i]

    def calculatePercentage(self):
        # Calculating total sum to find percentage
        categories_total = dict()
        biotypes_total = dict()
        for i in self.categories.keys():
            for j in self.categories[i]:
                categories_total[i] = categories_total.get(
                    i, 0) + self.categories[i][j]
        for i in self.biotypes.keys():
            for j in self.biotypes[i]:
                biotypes_total[i] = biotypes_total.get(
                    i, 0) + self.biotypes[i][j]

        # Divinding by total subiotypesm
        for i in self.categories.keys():
            for j in self.categories[i]:
                self.categories[i][j] = self.categories[i][j] / \
                    categories_total[i] * 100
        for i in self.biotypes.keys():
            for j in self.biotypes[i]:
                self.biotypes[i][j] = self.biotypes[i][j] / \
                    biotypes_total[i] * 100

    def print_alfa_charts(self):

        #cats = ['exon', 'exon']
        cats = []
        for i in self.categories.keys():
            for j in self.categories[i]:
                cats.append(str(j))
        # print(cats)

        config = {
            # Building the plot
            'id': '<random string>',                # HTML ID used for plot
            'cpswitch': False,                       # Show the 'Counts / Percentages' switch?
            # Initial display with 'Counts' specified? False for percentages.
            'cpswitch_c_active': True,
            'cpswitch_counts_label': 'Counts',      # Label for 'Counts' button
            'cpswitch_percent_label': 'Percentages',  # Label for 'Percentages' button
            'logswitch': False,                     # Show the 'Log10' switch?
            'logswitch_active': False,              # Initial display with 'Log10' active?
            'logswitch_label': 'Log10',             # Label for 'Log10' button
            # Hide categories where data for all samples is 0
            'hide_zero_cats': True,
            # Customising the plot
            # Plot title - should be in format "Module Name: Plot Title"
            'title': None,
            'xlab': "Sample",                           # X axis label
            # Y axis label
            'ylab': "Propotion of Reads(%)",
            'ymax': 100,                           # Max y limit
            'ymin': 0,                           # Min y limit
            # Maximum value for automatic axis limit (good for percentages)
            'yCeiling': None,
            'yFloor': None,                         # Minimum value for automatic axis limit
            'yMinRange': None,                      # Minimum range for axis
            'yDecimals': True,                      # Set to false to only show integer labels
            # Format string for x axis labels. Defaults to {value}
            'ylab_format': None,
            # Set to None to have category bars side by side
            'stacking': 'normal',
            'use_legend': True,                     # Show / hide the legend
            # Javascript function to be called when a point is clicked
            'click_func': None,
            'cursor': None,                         # CSS mouse cursor type.
            # Number of decimal places to use in the tooltip number
            'tt_decimals': 0,
            'tt_suffix': '',                        # Suffix to add after tooltip number
            # Show the percentages of each count in the tooltip
            'tt_percentages': False,
        }

        self.add_section(
            name='Categories',
            anchor='categories',
            plot=bargraph.plot(self.categories)
        )

        self.add_section(
            name='BioTypes',
            anchor='biotypes',
            plot=bargraph.plot(self.biotypes)
        )

        return
