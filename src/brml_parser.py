# Bruker .BRML format parser from PyXRD
#                                      PyXRD                                     
#     A python implementation of the matrix algorithm developed for the X-ray    
#         diffraction analysis of disordered lamellar structures         
#                      Copyright (c) 2013-2014, Mathijs Dumon                       
#        This software is licensed under a BSD-2 Clause ("FreeBSD") License,                              

import os
from zipfile import ZipFile
import xml.etree.ElementTree as ET
import numpy as np
from xml.etree import cElementTree as ElementTree

class BrmlParser:
    """
    Bruker .BRML format parser
    """

    def __init__(self, file_path):
        self.file_path = file_path
        self.data = None
        self.header_dict = None

        
    def _get_file(self):
        """
        Returns a tuple: (filename, zipfile-object)
        """
        return self.file_path, ZipFile(self.file_path, "r")

    def _get_raw_data_files(self, f, folder):
        """
        Processes DataContainer.xml and returns a list of xml raw data
        filepaths and the sample name
        """
        with f.open(f"{folder}/DataContainer.xml", "r") as contf:
            data = contf.read()

        root = ET.fromstring(data)
        sample_name = root.find("./MeasurementInfo").get("SampleName")

        raw_data_files = []
        for child in root.find("./RawDataReferenceList"):
            raw_data_files.append(child.text)

        return raw_data_files, sample_name

    def _get_header_dict(self, f, folder):
        header_d = {}

        with f.open(f"{folder}/MeasurementContainer.xml", "r") as contf:
            data = contf.read()

        root = ET.fromstring(data)

        radius_path = (
            "./HardwareLogicExt/Instrument/BeamPathContainers"
            + "/BeamPathContainerAbc[@VisibleName='PrimaryTrack']/%s"
        )

        tube_path = (
            "./HardwareLogicExt/Instrument/BeamPathContainers"
            + "/BeamPathContainerAbc[@VisibleName='PrimaryTrack']/BankPositions/"
            + "BankPosition/MountedComponent/MountedTube/%s"
        )

        soller1_path = (
            "./HardwareLogicExt/Instrument/BeamPathContainers"
            + "/BeamPathContainerAbc[@VisibleName='PrimaryTrack']/BankPositions/"
            + "BankPosition/MountedComponent[@VisibleName='SollerMount']/%s"
        )

        soller2_path = (
            "./HardwareLogicExt/Instrument/BeamPathContainers"
            + "/BeamPathContainerAbc[@VisibleName='SecondaryTrack']/BankPositions/"
            + "BankPosition/MountedComponent[@VisibleName='SollerMount']/%s"
        )

        divergence_path = (
            "./HardwareLogicExt/Instrument/BeamPathContainers"
            + "/BeamPathContainerAbc[@VisibleName='PrimaryTrack']/BankPositions/"
            + "BankPosition/MountedComponent/Restrictions[@FieldName='OpeningDegree']/%s"
        )

        header_d.update(
            alpha1=float(root.find(tube_path % "WaveLengthAlpha1").get("Value")),
            alpha2=float(root.find(tube_path % "WaveLengthAlpha2").get("Value")),
            alpha_average=float(root.find(tube_path % "WaveLengthAverage").get("Value")),
            beta=float(root.find(tube_path % "WaveLengthBeta").get("Value")),
            alpha_factor=root.find(tube_path % "WaveLengthRatio").get("Value"),
            target_type=root.find(tube_path % "TubeMaterial").get("Value"),
            soller1=root.find(soller1_path % "Deflection").get("Value"),
            soller2=root.find(soller2_path % "Deflection").get("Value"),
            radius=float(root.find(radius_path % "Radius").get("Value", 0)),
        )

        return header_d

    def parse(self):
        filename, f = self._get_file()

        folder = "Experiment0"
        # folder = os.path.splitext(os.path.basename(filename))[0]
        raw_data_files, sample_name = self._get_raw_data_files(f, folder)
        header_dict = self._get_header_dict(f, folder)
        zipinfos = f.infolist()

        data = dict()
        for raw_file in raw_data_files:
            two_theta = []
            intensity = []
            with f.open(f"{raw_file}", "r") as rawf:
                data_t = rawf.read()
            root = ET.fromstring(data_t)
            for route in root.findall("./DataRoutes/DataRoute"):
                datums = route.findall("Datum")
                data_t = []
                enabled_datum_index = None
                twotheta_datum_index = None
                intensity_datum_index = None
                steptime_datum_index = None
                for dataview in route.findall("./DataViews/RawDataView"):
                    index = int(dataview.get("Start", 0))
                    name = dataview.get("LogicName") or "Undefined"
                    xsi_type = (
                        dataview.get("{http://www.w3.org/2001/XMLSchema-instance}type")
                        or "Undefined"
                    )
                    if name == "MeasuredTime":
                        steptime_datum_index = index
                    elif name == "AbsorptionFactor":
                        enabled_datum_index = index
                    elif name == "Undefined" and xsi_type == "VaryingRawDataView":
                        for i, definition in enumerate(
                            dataview.findall("./Varying/FieldDefinitions")
                        ):
                            if definition.get("TwoTheta"):
                                index += i
                                break
                        twotheta_datum_index = index
                    elif name == "Undefined" and xsi_type == "RecordedRawDataView":
                        logic_name = dataview.find("Recording").get("LogicName")
                        if logic_name == "Counter0D" or logic_name == "ScanCounter":
                            intensity_datum_index = index
                        elif logic_name == "modeActualHum":
                            relative_humidity_index = index
                            relative_humidity_data = []
                        elif logic_name == "modeActualTemp":
                            temperature_index = index
                            temperature_data = []

            for subscan in route.findall("./SubScans/SubScanInfo"):
                # Get the steps, where to start and the planned
                # time per step (measuredTimePerStep deviates
                # if the recording was interrupted):
                steps = int(subscan.get("MeasuredSteps"))
                start = int(subscan.get("StartStepNo"))
                steptime = float(subscan.get("PlannedTimePerStep"))

                for datum in datums[start : start + steps]:
                    values = datum.text.split(",")
                    # Fetch values from the list:
                    datum_steptime = float(values[steptime_datum_index])
                    intensity = float(values[intensity_datum_index])
                    intensity /= float(steptime * datum_steptime)
                    twotheta = float(values[twotheta_datum_index])

                    # Append point and increase count:
                    data_t.append([twotheta, intensity])
            base = os.path.basename(raw_file).split(".")[0]
            data[base] = np.array(data_t)
        self.header_dict = header_dict
        self.data = data

    def get_data(self):
        if self.data is None:
            self.parse()
        return self.data, self.header_dict

    def save_to_csv(self, output_path):
        if self.data is None:
            self.parse()

        import numpy as np

        np.savetxt(
            output_path,
            np.column_stack((self.data["two_theta"], self.data["intensity"])),
            delimiter=",",
            header="2Theta,Intensity",
            comments="",
        )
