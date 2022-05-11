import argparse
from argparse import RawDescriptionHelpFormatter
from numpy import log10, floor, ceil, sqrt, arange
from tqdm import tqdm
import logging
import operator

class CalcSuppr():
    def __init__(self, freq, eff_suppr, field, fast_mode=False, precision=0.001):
        self.freq = freq * 10**9
        self.eff_suppr = eff_suppr
        self.field = field
        self.wave_length = (299792458/self.freq)*100
        self.max_line_length = floor((self.wave_length/2)*10000)/10000
        self.min_interval = ceil((self.wave_length/10)*10000)/10000
        self.area_interval = floor((self.wave_length/2)*10000)/10000
        self.fast_mode = fast_mode
        self.precision = precision

    def roundUp(self, number):
        return ceil(number*10000)/10000

    def roundDown(self, number):
        return floor(number*10000)/10000

    def calcSuppr(self, line_length, height, interval_x):
        if interval_x <= self.area_interval:
            close_apert = 6
            if interval_x*sqrt(3) + line_length/2 <= self.area_interval:
                close_apert += 6
                if 2*interval_x + height <= self.area_interval:
                    close_apert += 6
        else:
            close_apert = 0

        s = 20*log10(self.wave_length/(2*line_length)) - 20*log10(sqrt(close_apert+1))
        return s

    def numAperture(self, line_length, interval):
        height = self.roundUp((line_length/2)*sqrt(3))

        if self.fast_mode == True:
            interval_x = self.area_interval
            interval_y = self.roundUp(interval_x*sqrt(3)/2 - line_length/4)
        else:
            interval_x = interval
            interval_y = self.roundUp(interval*sqrt(3)/2 - line_length/4)
        
        temp_x_field = self.field[0] - 2*self.min_interval
        temp_y_field = self.field[1] - 2*self.min_interval

        num_rows = int(floor((temp_y_field+interval_y)/(line_length+interval_y)))

        counter=0
        fig_set = []
        for row in range(0, num_rows):
            if row%2 == 0:
                temp = int(floor((temp_x_field+interval_x)/(height+interval_x)))
                counter+= temp
                fig_set.append(temp)
            else:
                temp = int(floor((temp_x_field-height/2-interval_x/2+interval_x)/(height+interval_x)))
                counter += temp
                fig_set.append(temp)

        field_1fig = 3/2 * (line_length/2)**2 * sqrt(3)

        if self.fast_mode == True:
            s = 20*log10(self.wave_length/(2*line_length))
        else:
            s = self.calcSuppr(line_length, height, interval_x)


        return {"Field":counter * field_1fig, "Aperatures":counter, "Height":height,
        "Interval X":interval_x, "Interval Y":interval_y, "Suppression":s, "Set":fig_set}

    def getResults(self):
        logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
        results = [{"Line length":0, "Interval X":0, "Field":0}]
        line_length=0
        best_set = {"Field":0, "Aperatures":0}
        logging.info("Wave length: "+str(self.wave_length))
        logging.info(" Lambda/2: "+str(self.area_interval))
        logging.info("Lambda/10: "+str(self.min_interval))

        if self.fast_mode == True:
            for x in tqdm(arange(self.max_line_length/2, self.max_line_length-0.1, self.precision)):
                line_length=self.max_line_length-x
                aperatures_field = self.numAperture(line_length, 0)
                if aperatures_field["Suppression"] >= self.eff_suppr and aperatures_field["Suppression"] <= self.eff_suppr+0.3:
                        
                    if aperatures_field["Aperatures"] > best_set["Aperatures"] and \
                    aperatures_field["Field"] >= 0.01*best_set["Field"]:
                            
                        if line_length == results[-1]["Line length"] and \
                        aperatures_field["Field"] == results[-1]["Field"]:
                            continue

                        if aperatures_field["Aperatures"] > best_set["Aperatures"] and aperatures_field["Field"] > best_set["Field"]:
                            best_set = {"Field":aperatures_field["Aperatures"], "Aperatures":aperatures_field["Field"]}

                        temp_results = {
                            "Line length":line_length, "Field":aperatures_field["Field"], "Aperatures":aperatures_field["Aperatures"],
                            "Interval X": aperatures_field["Interval X"], "Interval Y": aperatures_field["Interval Y"], 
                            "Suppr":aperatures_field["Suppression"], "Fig set order per row": aperatures_field["Set"]
                        }
                        results.append(temp_results)
                        logging.info("\n")
                        logging.info(temp_results)
            
            results.pop(0)
            return sorted(results, key=operator.itemgetter("Field", "Aperatures"))       
            
        else:
            if self.field[0] > self.field[1]:
                max_interval = self.field[1]/3
            else:
                max_interval = self.field[0]/3

            for x in tqdm(arange(self.max_line_length/2, self.max_line_length-0.1, self.precision)):
                line_length=self.max_line_length-x

                for interval in arange(self.min_interval, max_interval, self.precision):
                    if self.area_interval > 2*interval+self.roundUp((line_length/2)*sqrt(3)):
                        continue
                    
                    aperatures_field = self.numAperture(line_length, interval)
                    if aperatures_field["Suppression"] >= self.eff_suppr and aperatures_field["Suppression"] <= self.eff_suppr+0.3:
                        
                        if aperatures_field["Aperatures"] > best_set["Aperatures"] and \
                        aperatures_field["Field"] >= 0.01*best_set["Field"]:
                            
                            if line_length == results[-1]["Line length"] and \
                            aperatures_field["Field"] == results[-1]["Field"]:
                                continue

                            if aperatures_field["Aperatures"] > best_set["Aperatures"] and aperatures_field["Field"] > best_set["Field"]:
                                best_set = {"Field":aperatures_field["Aperatures"], "Aperatures":aperatures_field["Field"]}

                            temp_results = {
                                "Line length":line_length, "Field":aperatures_field["Field"], "Aperatures":aperatures_field["Aperatures"],
                                "Height": aperatures_field["Height"], "Interval X": aperatures_field["Interval X"], "Interval Y": aperatures_field["Interval Y"], 
                                "Suppr":aperatures_field["Suppression"], "Fig set order per row": aperatures_field["Set"]
                            }
                            results.append(temp_results)
                            logging.info("\n")
                            logging.info(temp_results)
            
            results.pop(0)
            return sorted(results, key=operator.itemgetter("Field", "Aperatures"))    
        

parser = argparse.ArgumentParser(usage='use "python %(prog)s --help" for more information',formatter_class=RawDescriptionHelpFormatter, description='Created by Patryk Duniak 253003')
parser.add_argument('-f', '--frequency', type=float, required=True, help="Input int/float of frequenzy [GHz]")
parser.add_argument('-s', '--shielding_eff', type=float, required=True, help="Input int/float of desired shielding efectiveness [dB]")
parser.add_argument('-x', '--xaxis', type=float, required=True, help="Input int/float of field in x-axis [cm]")
parser.add_argument('-y', '--yaxis', type=float, required=True, help="Input int/float of field in y-axis [cm]")
parser.add_argument('-p', '--precision', type=float, default=0.01, help="Precision [cm] in iterations, may have an effect of speed. Default is 0.01" )
parser.add_argument('--fast', action='store_true', help="Fixed interval length")
args = parser.parse_args()

if args.fast == True:
    test = CalcSuppr(args.frequency, args.shielding_eff, [args.xaxis, args.yaxis], True, args.precision)
else:
    test = CalcSuppr(args.frequency, args.shielding_eff, [args.xaxis, args.yaxis], False, args.precision)

results = test.getResults()
for element in results:
    print("Wymiar Liniowy[cm]: "+str(element["Line length"]))
    print("Pole        [cm^2]: "+str(element["Field"]))
    print("Otwory            : "+str(element["Aperatures"]))
    print("Wysokosc fig. [cm]: "+str(element["Height"]))
    print("Odstep w osi x[cm]: "+str(element["Interval X"]))
    print("Odstep w osi y[cm]: "+str(element["Interval Y"]))
    print("Skutecznosc ekranowania     [dB]: "+str(element["Suppr"]))
    print("Ilosc otworow w danym rzedzie :"+str(element["Fig set order per row"]))
    print("\n")