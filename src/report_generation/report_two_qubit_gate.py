######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: report_two_qubit_gate.py

# Author(s): Henrik Zander

# Date created: 19 December 2021

# Copyright 2021-2022, Henrik Zander, All rights reserved.

######################################################################################################################################################################

# Global Python imports
from fpdf import FPDF
from pathlib import Path
import os

######################################################################################################################################################################

# Variable initializations
debug = False
WIDTH = 210
HEIGHT = 297

# Path initializations
graphDirectory = Path.cwd() / "temp"
symbolsDirectory = Path.cwd() / "report_generation" / "Symbols"
brandDirectory = Path.cwd() / "report_generation" / "Branding"
equationsDirectory = Path.cwd() / "report_generation" / "Equations"

######################################################################################################################################################################
# Class definitions

# Redefinition of the "footer"-function so that page number can be added.
class PDF(FPDF):
    def footer(self):
        # Set position of the footer
        self.set_y(-15)
        # Set font of footer
        self.set_font('helvetica', 'I', 10)
        # Add page number
        self.cell(0, 10, f'Page {self.page_no()}/{{nb}}', align='C', border=debug)


######################################################################################################################################################################
# Helper functions

def numOfAvailableGraphs():
    return len([name for name in os.listdir(graphDirectory)])


def insertReportTitle(pdf, title, fontSize):
    pdf.set_font('helvetica', 'B', fontSize)
    pdf.ln(50)
    pdf.cell(0, 15, title, border=debug, ln=True)


def insertTitle(pdf, title, fontSize):
    pdf.set_font('helvetica', 'UB', fontSize)
    pdf.cell(0, 10, title, border=debug, ln=True)


def insertGeneralParameters(pdf, result):
    # Change back the body text font
    pdf.set_font('helvetica', '', 12)

    # Create the body text for the 'General Parameters'-section
    body_txt = f'Chip Name: e.g. the 5 qubit chip, unless we have codenames\nTargeted Qubits: e.g. Q1, Q4\nDate and Time: {result["creationTime"]}\n\nGate: {result["gateType"]}\nFidelity: {round(result["gateFidelity"], 5)}\nModulation Time: {round(result["modulationTime"], 5)} ns\nEnergy Levels per Qubit: {result["nOptimizationLvls"]}'

    # Insert the body text for the 'General Parameters'-section
    pdf.multi_cell(0, 5, body_txt, border=debug)
    pdf.ln(10)


def insertSignalParameters(pdf, result):
    # Create the body text for the 'Signal Parameters'-section
    if result["signalType"] == 'cos':

        # Make the text underlined
        pdf.set_font('helvetica', 'U', 12)

        # Insert picture of the signal shape
        pdf.cell(0, 5, "Signal Shape:", ln=True, border=debug)
        pdf.image(str((equationsDirectory / "cos-signal.png").resolve()), w = 60)
        pdf.ln(5)

        """# Insert picture of tunable bus frequency
        pdf.cell(0, 5, "Approximate Frequency of the Tunable Bus:", ln=True, border=debug)
        pdf.image(str((equationsDirectory / "tunable-bus-frequency.png").resolve()), w = 60)
        pdf.ln(4) """

        # Insert picture of amplitude modulation
        pdf.cell(0, 5, "Amplitude Modulation:", ln=True, border=debug)
        pdf.image(str((equationsDirectory / "delta_modulation.png").resolve()), x=11, w = 32)
        pdf.ln(5)

        # Insert title for the 'Optimal Signal'-subsection
        pdf.cell(0, 5, "Parameters of Optimal Signal:", ln=True, border=debug)
        pdf.ln(2)

        # Create the body text
        body_txt = f'   : {round(result["theta"], 5)} [      ]\n   : {round(result["delta"], 5)} [      ]\n     : {round(result["omegaPhi"], 5)} GHz\nTotal Modulation Time of       : {round(result["modulationTime"], 5)} ns\nRise Time of Modulation (0 to 100%): {result["riseTime"]} ns'

        # Insert the latex pictures of the variables
        pdf.image(str((symbolsDirectory / "Theta.png").resolve()), x = 11, y = 197.5, w = 4)
        pdf.image(str((symbolsDirectory / "Phi_0.png").resolve()), x = 36.5, y = 197.5, w = 6)
        pdf.image(str((symbolsDirectory / "delta_0.png").resolve()), x = 11, y = 202.5, w = 4)
        pdf.image(str((symbolsDirectory / "Phi_0.png").resolve()), x = 35, y = 202.5, w = 6)
        pdf.image(str((symbolsDirectory / "omega_Phi.png").resolve()), x = 11, y = 208, w = 6)
        pdf.image(str((symbolsDirectory / "delta_t.png").resolve()), x = 59, y = 212, w = 7)

    elif result["signalType"] == 'arccos':
        pass # ToDO: Add code for the arccos signal

    # Change back the body text font
    pdf.set_font('helvetica', '', 12)

    # Insert the body text for the 'Signal Parameters'-section
    pdf.multi_cell(0, 5, body_txt, border=debug)
    pdf.ln(10)


def insertCircuitParameters(pdf, result):
    # Make the text underlined
    pdf.set_font('helvetica', 'U', 12)

    # Insert the "qubit frequencies"-title.
    pdf.cell(0, 5, "Frequencies:", ln=True, border=debug)

    # Change back the body text font
    pdf.set_font('helvetica', '', 12)

    # Create the "frequencies"-text.
    freq_text = f'\tQ1: {result["frequencies"][0]} GHz\n\tQ2: {result["frequencies"][1]} GHz\n\tCoupler: {result["frequencies"][2]} GHz'
    
    # Insert the "frequencies"-text.
    pdf.multi_cell(0, 5, freq_text, border=debug)
    pdf.ln(5)

    # Add a new page so that the page-break is better looking.
    pdf.add_page()

    # Make the text underlined
    pdf.set_font('helvetica', 'U', 12)

    # Insert the "qubit anharmonicities"-title.
    pdf.cell(0, 5, "Anharmonicities:", ln=True, border=debug)

    # Change back the body text font
    pdf.set_font('helvetica', '', 12)

    # Create the "anharmonicities"-text.
    anharm_text = f'\tQ1: {result["anharmonicities"][0]} GHz\n\tQ2: {result["anharmonicities"][1]} GHz\n\tCoupler: {result["anharmonicities"][2]} GHz'
    
    # Insert the "anharmonicities"-text.
    pdf.multi_cell(0, 5, anharm_text, border=debug)
    pdf.ln(5)

    # Make the text underlined
    pdf.set_font('helvetica', 'U', 12)

    # Insert the "coupling strengths"-title.
    pdf.cell(0, 5, "Coupling Strengths:", ln=True, border=debug)

    # Change back the body text font
    pdf.set_font('helvetica', '', 12)

    # Create the "coupling strengths"-text.
    coupling_strength_text = f'\tQ1 to Coupler: {result["couplings"][0]} GHz\n\tQ2 to Coupler: {result["couplings"][1]} GHz'
    
    # Insert the "coupling strengths"-text.
    pdf.multi_cell(0, 5, coupling_strength_text, border=debug)
    pdf.ln(10)


def insertGraphs(pdf):
    pdf.image(str((graphDirectory / "fig1.jpeg").resolve()), w=(WIDTH-30)/2, x=10, y=80)
    pdf.image(str((graphDirectory / "fig2.jpeg").resolve()), w=(WIDTH-30)/2, x=(WIDTH/2)+5, y=80)
    pdf.image(str((graphDirectory / "fig3.jpeg").resolve()), w=(WIDTH-30)/2, x=10, y=180)
    pdf.image(str((graphDirectory / "fig4.jpeg").resolve()), w=(WIDTH-30)/2, x=(WIDTH/2)+5, y=180)
    
    pdf.add_page()

    pdf.image(str((graphDirectory / "fig5.jpeg").resolve()), w=(WIDTH-30)/2, x=10, y=20)
    pdf.image(str((graphDirectory / "fig6.jpeg").resolve()), w=(WIDTH-30)/2, x=(WIDTH/2)+5, y=20)
    pdf.image(str((graphDirectory / "fig7.jpeg").resolve()), w=(WIDTH-30)/2, x=10, y=110)
    pdf.image(str((graphDirectory / "fig8.jpeg").resolve()), w=(WIDTH-30)/2, x=(WIDTH/2)+5, y=110)
    pdf.image(str((graphDirectory / "fig9.jpeg").resolve()), w=(WIDTH-30)/2, x=10, y=200)


######################################################################################################################################################################
# Optimization Result Report Generation

def generateReport(result):
    # Create the PDF object
    pdf = PDF('P', 'mm', 'A4')

    # Set nb to an alias for the total number of pages
    pdf.alias_nb_pages()
    # Change page break settings
    pdf.set_auto_page_break(auto=True, margin=15)
    # Add a page to the pdf object
    pdf.add_page()

    # Insert the letter-head picture
    pdf.image( str((brandDirectory / "letterHead.jpg").resolve()), w=WIDTH, x=0, y=0)

    # Insert the logo
    pdf.image( str((brandDirectory / "GATESIDE_BLACK_BLUE_BACKGROUND.jpg").resolve()), w=60, x=WIDTH-65, y=5)
    
    # Insert the main title of the report
    insertReportTitle(pdf, 'Optimization Report', 24)
    pdf.ln(10)

    # Insert the title for the 'General Parameters'-section
    insertTitle(pdf, 'General Parameters', 18)
    
    # Insert the "General Parameters"-section of the report
    insertGeneralParameters(pdf, result)

    # Insert the title for the 'Signal Parameters'-section
    insertTitle(pdf, 'Signal Parameters', 18)

    # Insert the "Signal Parameters"-section of the report
    insertSignalParameters(pdf, result)

    # Insert the title for the 'Signal Parameters'-section
    insertTitle(pdf, 'Circuit Parameters', 18)

    # Insert the "Signal Parameters"-section of the report
    insertCircuitParameters(pdf, result)
    
    # Insert the title for the 'Graphs'-section
    insertTitle(pdf, 'Graphs', 18)

    # Insert the graphs
    insertGraphs(pdf)

    # Output the PDF-file
    pdf.output('report.pdf')


######################################################################################################################################################################