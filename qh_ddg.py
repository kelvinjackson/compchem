#!/usr/bin/python

# THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# Comments and/or additions are welcome (send e-mail to:
# robert.paton@chem.ox.ac.uk

#######################################################################
#                               qh_ddg.py                             #
#  A program to recompute the vibrational entropy from a standard     #
#  output file as produced by a frequency calculation in Gaussian 09  #
#  The harmonic approximation is used for frequencies above a certain #
#  cut-off frequency, while the free-rotor approximation is applied   #
#  below this value. A damping function feathers between these two    #
#  expressions for the Svib. rather than a hard cut-off. This avoids  #
#  values of Svib. that tend to infinite values as the frequency      #
#  tends to zero. This approach was described in: Grimme, S. Chem.    #
#  Eur. J. 2012, 18, 9955                                             #
#
#  The free energy is then reevaluated at a specified temperature,    #
#  with a specified frequency cut-off. An optional vibrational        #
#  scaling factor and concentration may be specified. The latter term #
#  is factored into calculation of the translational entropy. Default #
#  values for these terms is 1.0 and 1 mol/l                          #
#######################################################################
#######  Written by:  Rob Paton #######################################
#######  Last modified:  Mar 20, 2013 #################################
#######################################################################

import sys, math

from GoodVibes import *

if __name__ == "__main__":
	
	# Takes arguments: cutoff_freq g09_output_files
	files = []
	FREQ_CUTOFF = "none"; temperature = "none"; conc = "none"; freq_scale_factor = "none"; solv = "none"
	if len(sys.argv) > 1:
		for i in range(1,len(sys.argv)):
			if sys.argv[i] == "-f": FREQ_CUTOFF = float(sys.argv[i+1])
			elif sys.argv[i] == "-t": temperature = float(sys.argv[i+1])
			elif sys.argv[i] == "-c": conc = float(sys.argv[i+1])
			elif sys.argv[i] == "-v": freq_scale_factor = float(sys.argv[i+1])
			elif sys.argv[i] == "-s": solv = (sys.argv[i+1])
			
			else:
				if len(sys.argv[i].split(".")) > 1:
					if sys.argv[i].split(".")[1] == "out" or sys.argv[i].split(".")[1] == "log": files.append(sys.argv[i])
	
		if FREQ_CUTOFF != "none": print "   Frequency cut-off value =", FREQ_CUTOFF, "wavenumbers"
		else: print "   Frequency cut-off value not defined!"; sys.exit()
		if temperature != "none": print "   Temperature =", temperature, "Kelvin",
		else: print "   Temperature (default) = 298.15K",; temperature = 298.15
		if conc != "none": print "   Concn =", conc, "mol/l",
		else: print "   Concn (default) = 1 mol/l",; conc = 1.0
		if freq_scale_factor != "none": print "   Frequency scale factor =", freq_scale_factor,
		else: print "   Frequency scale factor (default) = 1.0"; freq_scale_factor = 1.0
		
		print solv
		freespace = get_free_space(solv)
		if freespace != 1000.0: print "   Solvent =", solv+": % free volume (Strans)","%.1f" % (freespace/10.0)
		else: print "   Solvent (default) = undefined: % free volume (Strans) = 100.0"; solv = "unk"
	
	else:
		print "\nWrong number of arguments used. Correct format: qh_ddg.py -f cutoff_freq (-t temp) (-c concn) (-s scalefactor) g09_output_files\n"
		sys.exit()

	print "\n  ",
	print "Structure".ljust(40), "Energy".rjust(10), "ZPE".rjust(10), "Enthalpy".rjust(10), "T.S".rjust(10), "T.qh-S".rjust(10), "G(T)".rjust(10), "qh-G(T)".rjust(10)
	ddg = []
	for file in files:
		bbe = calc_bbe(file, FREQ_CUTOFF, temperature, conc, freq_scale_factor, solv)
		ddg.append(bbe)
		print "o ",
		print (file.split(".")[0]).ljust(40),
		if not hasattr(bbe,"gibbs_free_energy"): print "Warning! Job did not finish normally!"
		if hasattr(bbe, "scf_energy"): print "%.6f" % bbe.scf_energy,
		else: print "N/A",
		if hasattr(bbe, "zero_point_corr"): print "   %.6f" % (bbe.zpe),
		else: print "N/A",
		if hasattr(bbe, "enthalpy"): print "   %.6f" % (bbe.enthalpy),
		else: print "N/A",
		if hasattr(bbe, "entropy"): print "   %.6f" % (temperature * bbe.entropy),
                else: print "N/A",
		if hasattr(bbe, "qh_entropy"): print "   %.6f" % (temperature * bbe.qh_entropy),
                else: print "N/A",
		if hasattr(bbe, "gibbs_free_energy"): print "   %.6f" % (bbe.gibbs_free_energy),
		else: print "N/A",
		if hasattr(bbe, "qh_gibbs_free_energy"): print "   %.6f" % (bbe.qh_gibbs_free_energy)
                else: print "N/A"
	print ""
	dde = autokcal * (ddg[0].scf_energy - ddg[1].scf_energy)
	ddzpe = autokcal * (ddg[0].zpe - ddg[1].zpe)
	ddh = autokcal * (ddg[0].enthalpy - ddg[1].enthalpy)
	dds = autokcal * temperature * (ddg[0].entropy - ddg[1].entropy)
	ddqhs = autokcal * temperature * (ddg[0].qh_entropy - ddg[1].qh_entropy)
	ddgibbs = autokcal * (ddg[0].gibbs_free_energy - ddg[1].gibbs_free_energy)
	ddqhgibbs = autokcal * (ddg[0].qh_gibbs_free_energy - ddg[1].qh_gibbs_free_energy)
	print "".ljust(40), "%.6f" % dde,  "   %.6f" % (ddzpe), "   %.6f" % (ddh), "   %.6f" % (dds), "   %.6f" % (ddqhs), "   %.6f" % (ddgibbs), "   %.6f" % (ddqhgibbs),
