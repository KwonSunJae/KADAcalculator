from email.contentmanager import raw_data_manager
from combus import CommunicationBus
import convert
from flight_conditions import FlightConditions
import paths
import sys
import numpy as np
from scipy.interpolate import interp1d
import subprocess as sp
import os
import shutil

import matplotlib.pyplot as plt


def get_airfoil_thickness(filename):
    fid = open('./ram-analysis/'+filename, 'rt')
    lines = fid.readlines()
    fid.close()
    npts = len(lines)
    pts = np.zeros([npts, 2], float)
    for i,line in enumerate(lines):
        val = line.split()
        pts[i] = np.array([float(val) for val in line.split()])
    return np.max(pts[:,1]) - np.min(pts[:,1])


class IntegratedAnalysis(object):
    def __init__(self) -> None:
        self.cb = CommunicationBus('combus_ram.csv')
        self.read_input()
        self.analyze_geometry()
    
    def read_input(self):
        f = open( sys.argv[1], 'r')
        lines = f.readlines()

        inpData = []
        for line in lines : 
            item = line.split(" ")
            inpData.append( item )

        itr = len(inpData)
        for i in range(0, itr):
            if 'fus_apex_FS' in inpData[i]:
                self.cb.fus_apex_FS = float(inpData[i][2])
            if 'fus_apex_WL' in inpData[i]:
                self.cb.fus_apex_WL = float(inpData[i][2])
            if 'fus_diameter' in inpData[i]:
                self.cb.fus_diameter = float(inpData[i][2])
            if 'fus_length' in inpData[i]:
                self.cb.fus_length = float(inpData[i][2])
            if 'fus_nose_to_length_ratio' in inpData[i]:
                self.cb.fus_nose_to_length_ratio = float(inpData[i][2])
            if 'fus_tail_to_length_ratio' in inpData[i]: 
                self.cb.fus_tail_to_length_ratio = float(inpData[i][2])

            if 'wing_area' in inpData[i]:
                self.cb.wing_area = float(inpData[i][2])
            if 'wing_aspect_ratio' in inpData[i]:
                self.cb.wing_aspect_ratio = float(inpData[i][2])
            if 'wing_taper_ratio' in inpData[i]:
                self.cb.wing_taper_ratio = float(inpData[i][2])
            if 'wing_leading_edge_sweep' in inpData[i]:
                self.cb.wing_leading_edge_sweep = float(inpData[i][2])
            if 'wing_apex_FS' in inpData[i]:
                self.cb.wing_apex_FS = float(inpData[i][2])
            if 'wing_apex_WL' in inpData[i]:
                self.cb.wing_apex_WL = float(inpData[i][2])
            if 'wing_incidence_angle' in inpData[i]:
                self.cb.wing_incidence_angle = float(inpData[i][2])
            if 'wing_airfoil_filename' in inpData[i]:
                self.cb.wing_airfoil_filename = str(inpData[i][2][:-1])
            
            if 'hstab_area' in inpData[i]:
                self.cb.hstab_area = float(inpData[i][2])
            if 'hstab_aspect_ratio' in inpData[i]:
                self.cb.hstab_aspect_ratio = float(inpData[i][2])
            if 'hstab_taper_ratio' in inpData[i]:
                self.cb.hstab_taper_ratio = float(inpData[i][2])
            if 'hstab_leading_edge_sweep' in inpData[i]:
                self.cb.hstab_leading_edge_sweep = float(inpData[i][2])
            if 'hstab_apex_FS' in inpData[i]:
                self.cb.hstab_apex_FS = float(inpData[i][2])
            if 'hstab_apex_WL' in inpData[i]:
                self.cb.hstab_apex_WL = float(inpData[i][2])
            if 'hstab_incidence_angle' in inpData[i]:
                self.cb.hstab_incidence_angle = float(inpData[i][2])
            if 'hstab_airfoil_filename' in inpData[i]:
                self.cb.hstab_airfoil_filename = str(inpData[i][2][:-1])
            if 'elevator_to_chord_ratio' in inpData[i]:
                self.cb.elevator_to_chord_ratio = float(inpData[i][2])
            
            if 'vstab_area' in inpData[i]:
                self.cb.vstab_area = float(inpData[i][2])
            if 'vstab_aspect_ratio' in inpData[i]:
                self.cb.vstab_aspect_ratio = float(inpData[i][2])
            if 'vstab_taper_ratio' in inpData[i]:
                self.cb.vstab_taper_ratio = float(inpData[i][2])
            if 'vstab_leading_edge_sweep' in inpData[i]:
                self.cb.vstab_leading_edge_sweep = float(inpData[i][2])
            if 'vstab_apex_FS' in inpData[i]:
                self.cb.vstab_apex_FS = float(inpData[i][2])
            if 'vstab_apex_WL' in inpData[i]:
                self.cb.vstab_apex_WL = float(inpData[i][2])
            if 'vstab_airfoil_filename' in inpData[i]:
                self.cb.vstab_airfoil_filename = str(inpData[i][2][:-1])
            
            if 'mass' in inpData[i]:
                self.cb.mass = float(inpData[i][2])
            if 'cg_percent_of_mac' in inpData[i]:
                self.cb.cg_percent_of_mac = float(inpData[i][2])
            
            if 'velocity' in inpData[i]:
                self.cb.velocity = float(inpData[i][2])
            if 'altitude' in inpData[i]:
                self.cb.altitude = float(inpData[i][2])
            if 'aoa_min' in inpData[i]:
                self.cb.aoa_min = float(inpData[i][2])
            if 'aoa_max' in inpData[i]:
                self.cb.aoa_max = float(inpData[i][2])
            if 'aoa_step' in inpData[i]:
                self.cb.aoa_step = float(inpData[i][2])
            
            if 'velocity_min' in inpData[i]:
                self.cb.velocity_min = float(inpData[i][2])
            if 'velocity_max' in inpData[i]:
                self.cb.velocity_max = float(inpData[i][2])
            if 'velocity_step' in inpData[i]:
                self.cb.velocity_step = float(inpData[i][2])

    def set_flight_condition(self):
        fc = FlightConditions(self.cb.velocity, self.cb.altitude)
        self.cb.Mach = fc.Mach
        self.cb.air_density = fc.atm.density
        self.cb.grav_accel = fc.g
        self.cb.dynamic_pressure = fc.dynamic_pressure
        self.cb.speed_of_sound = fc.atm.sound_speed

    def analyze_geometry(self):
        # airfoil
        self.cb.wing_thickness = \
            get_airfoil_thickness(self.cb.wing_airfoil_filename)
        self.cb.hstab_thickness = \
            get_airfoil_thickness(self.cb.hstab_airfoil_filename)
        self.cb.vstab_thickness = \
            get_airfoil_thickness(self.cb.vstab_airfoil_filename)
        # wing
        wtc = self.cb.wing_thickness
        htc = self.cb.hstab_thickness
        vtc = self.cb.vstab_thickness
        self.cb.wing_wetted_area = self.cb.wing_area *\
            (2+0.1843*wtc + 1.5268*wtc**2 - 0.8395*wtc**3)
        self.cb.hstab_wetted_area = self.cb.hstab_area *\
            (2+0.1843*htc + 1.5268*htc**2 - 0.8395*htc**3)
        self.cb.vstab_wetted_area = self.cb.vstab_area *\
            (2+0.1843*vtc + 1.5268*vtc**2 - 0.8395*vtc**3)
        
        self.cb.wing_span = np.sqrt(self.cb.wing_aspect_ratio\
            *self.cb.wing_area)
        self.cb.wing_root_chord = 2*self.cb.wing_span \
            / (self.cb.wing_aspect_ratio*(1+self.cb.wing_taper_ratio))
        self.cb.wing_tip_chord = self.cb.wing_root_chord * \
            self.cb.wing_taper_ratio

        self.cb.hstab_span = np.sqrt(self.cb.hstab_aspect_ratio\
            *self.cb.hstab_area)
        self.cb.hstab_root_chord = 2*self.cb.hstab_span \
            / (self.cb.hstab_aspect_ratio*(1+self.cb.hstab_taper_ratio))
        self.cb.hstab_tip_chord = self.cb.hstab_root_chord * \
            self.cb.hstab_taper_ratio

        self.cb.vstab_span = np.sqrt(self.cb.vstab_aspect_ratio\
            *self.cb.vstab_area)
        self.cb.vstab_root_chord = 2*self.cb.vstab_span \
            / (self.cb.vstab_aspect_ratio*(1+self.cb.vstab_taper_ratio))
        self.cb.vstab_tip_chord = self.cb.vstab_root_chord * \
            self.cb.vstab_taper_ratio
        
        wtr = self.cb.wing_taper_ratio
        wcr = self.cb.wing_root_chord
        self.cb.wing_mac = 2*wcr * (1 + wtr + wtr*wtr)/(3*(1 + wtr))
        self.cb.wing_mac_y_loc = self.cb.wing_span/6 * (1+2*wtr)/(1+wtr)
        self.cb.wing_tan_le = np.tan(np.deg2rad(
            self.cb.wing_leading_edge_sweep))        
        
        self.cb.aircraft_cgx = self.cb.wing_apex_FS \
            + self.cb.wing_mac_y_loc \
            *  self.cb.wing_tan_le + self.cb.cg_percent_of_mac \
            * self.cb.wing_mac
        self.cb.aircraft_cgz = self.cb.fus_apex_WL

        htr = self.cb.hstab_taper_ratio
        hcr = self.cb.hstab_root_chord
        self.cb.hstab_mac = 2*hcr * (1 + htr + htr*htr)/(3*(1 + htr))
        self.cb.hstab_tan_le = np.tan(np.deg2rad(
            self.cb.hstab_leading_edge_sweep))

        vtr = self.cb.vstab_taper_ratio
        vcr = self.cb.vstab_root_chord
        self.cb.vstab_mac = 2*vcr * (1 + vtr + vtr*vtr)/(3*(1 + vtr))
        self.cb.vstab_tan_le = np.tan(np.deg2rad(
            self.cb.vstab_leading_edge_sweep))

        # fuselage
        fus_cs_area = 0.25*np.pi*self.cb.fus_diameter**2.0
        fus_nose = self.cb.fus_nose_to_length_ratio * self.cb.fus_length
        fus_tail = self.cb.fus_tail_to_length_ratio * self.cb.fus_length
        self.cb.fus_wetted_area = (2.8*fus_nose + 2.5*fus_tail\
            + 4*(self.cb.fus_length - fus_nose - fus_tail)) \
            * np.sqrt(0.25*np.pi*fus_cs_area)
        self.cb.fus_diam_to_len_ratio = self.cb.fus_diameter \
            /self.cb.fus_length
        
        # fuselage profile
        fus_mid_section = self.cb.fus_length - fus_nose - fus_tail
        n = 21
        x = np.zeros(n)
        y = np.zeros(n)
        angle = np.flip(np.linspace(np.pi/2, np.pi, 7))
        x[:7] = np.array([fus_nose*np.cos(ang) for ang in angle]) \
             + fus_nose
        y[:7] = np.array([0.5*self.cb.fus_diameter*np.sin(ang) \
            for ang in angle])
        ptsx = np.array([fus_nose, fus_nose+fus_mid_section, 
            self.cb.fus_length])
        ptsy = np.array([self.cb.fus_diameter/2, self.cb.fus_diameter/2, 0])

        x[7:21] = np.linspace(ptsx[0], ptsx[2], 15, endpoint=True)[1:]
        f = interp1d(ptsx, ptsy, kind='linear')
        y[7:21] = f(x[7:21])

        x = np.hstack([np.flip(x), x[1:]])
        y = np.hstack([np.flip(y), -y[1:]])
        x += self.cb.fus_apex_FS
        y += self.cb.fus_apex_WL
        
        fid = open(paths.avl_fus_profile, 'wt')
        for _x, _y in zip(x,y):
            fid.write('%.6f  %.6f\n'%(_x, _y))
        fid.close()

    def calc_friction_drag(self):
        inp_file = os.path.abspath(paths.friction_inp)

        names = list()
        s_wet = list()
        l_ref = list()
        thickness = list()
        comp_type = list()

        names.append('WING')
        s_wet.append(convert.sqm_to_sqft(self.cb.wing_wetted_area))
        l_ref.append(convert.m_to_ft(self.cb.wing_mac))
        thickness.append(self.cb.wing_thickness)
        comp_type.append(0)

        names.append('HTAIL')
        s_wet.append(convert.sqm_to_sqft(self.cb.hstab_wetted_area))
        l_ref.append(convert.m_to_ft(self.cb.hstab_mac))
        thickness.append(self.cb.hstab_thickness)
        comp_type.append(0)

        names.append('VTAIL')
        s_wet.append(convert.sqm_to_sqft(self.cb.vstab_wetted_area))
        l_ref.append(convert.m_to_ft(self.cb.vstab_mac))
        thickness.append(self.cb.vstab_thickness)
        comp_type.append(0)

        names.append('BODY')
        s_wet.append(convert.sqm_to_sqft(self.cb.fus_wetted_area))
        l_ref.append(convert.m_to_ft(self.cb.fus_diameter))
        thickness.append(self.cb.fus_diam_to_len_ratio)
        comp_type.append(1)

        finp = open(inp_file, 'wt')

        finp.write('RAM\n')
        fmt = '{0:<10.4f}{1:<11.1f}{2:<10.1f}{3:<10.1f}\n'
        finp.write(fmt.format(convert.sqm_to_sqft(self.cb.wing_area), 
                   1, len(names), 0))

        fmt = '{0:<20s}{1:<10.4f}{2:<10.4f}{3:<10.6f}{4:<10.1f}{5:<10.4f}\n'

        for n,s,l,t,c in zip(names, s_wet, l_ref, thickness, comp_type):
            finp.write(fmt.format(n, s, l, t, c, 0))
        
        finp.write('{0:<10f}{1:<10f}\n'.format(self.cb.Mach, 
            convert.m_to_ft(self.cb.altitude)/1e3))
        finp.close()

        ps=sp.Popen(paths.friction, stdin=sp.PIPE, stdout=sp.PIPE, 
                    stderr=sp.PIPE)  
        cmd = 'default\n'
        ps.stdin.write(cmd.encode('utf-8'))
        cmd = paths.friction_inp + '\n'
        ps.stdin.write(cmd.encode('utf-8'))
        cmd = paths.friction_out + '\n'
        ps.stdin.write(cmd.encode('utf-8'))
        cmd = '\n'
        ps.stdin.write(cmd.encode('utf-8'))
        ps.communicate(('\n').encode('utf-8'))
        ps.kill()

        out_file = os.path.abspath(paths.friction_out)
        fid = open(out_file,'rt')
        lines = fid.readlines()
        fid.close()
        os.remove(paths.friction_inp)
        os.remove(paths.friction_out)
        self.cb.CD0 = float(lines[35].split()[-1])

    def create_avl_input_file(self):
        inp_filepath = os.path.abspath(paths.avl_inp)
        fid = open(inp_filepath,'wt')
        fid.write('RAM\n')
        fid.write('%.4f\n'%self.cb.Mach)
        fid.write('0  0  0\n')
        fid.write('%.4f  %.4f  %.4f\n'%(self.cb.wing_area, 
            self.cb.wing_mac, self.cb.wing_span))
        fid.write('%.4f  0.0  %.4f\n'%(self.cb.aircraft_cgx,
            self.cb.aircraft_cgz))
        fid.write('%.6f\n'%self.cb.CD0)
        
        fid.write('BODY\nfuselage\n')
        fid.write('28  1\n')
        fid.write('BFIL\n')
        fid.write('%s\n'%paths.avl_fus_profile)
        fid.write('SURFACE\n')
        fid.write('WING\n')
        fid.write('10  0  20  0\n')
        fid.write('YDUPLICATE\n0\n')
        # apex, chord, angle
        xw1 = self.cb.wing_apex_FS
        zw1 = self.cb.wing_apex_WL
        cw1 = self.cb.wing_root_chord
        iw = self.cb.wing_incidence_angle
        xw2 = xw1 + 0.5*self.cb.wing_span*self.cb.wing_tan_le
        yw2 = 0.5*self.cb.wing_span
        zw2 = zw1
        cw2 = self.cb.wing_tip_chord
        fid.write('SECTION\n')
        fid.write(f'{xw1}  0  {zw1}  {cw1}  {iw}\n')
        fid.write(f'AFIL\n{self.cb.wing_airfoil_filename}\n')
        fid.write('SECTION\n')
        fid.write(f'{xw2}  {yw2}  {zw2}  {cw2}  {iw}\n')
        fid.write(f'AFIL\n{self.cb.wing_airfoil_filename}\n')

        fid.write('SURFACE\n')
        fid.write('HSTAB\n')
        fid.write('6  0  10  0\n')
        fid.write('YDUPLICATE\n0\n')
        xh1 = self.cb.hstab_apex_FS
        zh1 = self.cb.hstab_apex_WL
        ch1 = self.cb.hstab_root_chord
        ih = self.cb.hstab_incidence_angle
        xh2 = xh1 + 0.5*self.cb.hstab_span*self.cb.hstab_tan_le
        yh2 = 0.5*self.cb.hstab_span
        zh2 = zh1
        ch2 = self.cb.hstab_tip_chord
        fid.write('SECTION\n')
        fid.write(f'{xh1}  0  {zh1}  {ch1}  {ih}\n')
        fid.write(f'AFIL\n{self.cb.hstab_airfoil_filename}\n')
        fid.write('CONTROL\n')
        fid.write('elevator  1.0  %.6f  0  0  0  +1\n'%\
            (1-self.cb.elevator_to_chord_ratio))
        fid.write('SECTION\n')
        fid.write(f'{xh2}  {yh2}  {zh2}  {ch2}  {ih}\n')
        fid.write(f'AFIL\n{self.cb.hstab_airfoil_filename}\n')
        fid.write('CONTROL\n')
        fid.write('elevator  1.0  %.6f  0  0  0  +1\n'%\
            (1-self.cb.elevator_to_chord_ratio))
        
        fid.write('SURFACE\n')
        fid.write('VSTAB\n')
        fid.write('6  0  10  0\n')
        xv1 = self.cb.vstab_apex_FS
        zv1 = self.cb.vstab_apex_WL
        cv1 = self.cb.vstab_root_chord
        xv2 = xv1 + self.cb.vstab_span*self.cb.vstab_tan_le
        zv2 = zv1 + self.cb.vstab_span
        cv2 = self.cb.vstab_tip_chord
        fid.write('SECTION\n')
        fid.write(f'{xv1}  0  {zv1}  {cv1}  0\n')
        fid.write(f'AFIL\n{self.cb.vstab_airfoil_filename}\n')
        fid.write('SECTION\n')
        fid.write(f'{xv2}  0  {zv2}  {cv2}  0\n')
        fid.write(f'AFIL\n{self.cb.vstab_airfoil_filename}\n')
        fid.close()

    def calc_avl_aero(self, aoa):
        ps=sp.Popen(paths.avl, stdin=sp.PIPE, stdout=sp.PIPE, 
                    stderr=sp.PIPE)
        cmd = 'load  %s\n'%paths.avl_inp  
        cmd += 'oper\na\na\n%.6f\n'%aoa
        cmd += 'm\nmn\n%.6f\n'%self.cb.Mach
        cmd += 'v\n%.6f\n'%self.cb.velocity
        cmd += 'd\n%.6f\n'%self.cb.air_density
        cmd += 'g\n%.6f\n'%self.cb.grav_accel
        cmd += 'm\n%.6f\n\n'%self.cb.mass
        cmd += 'x\nst\n%s\n\nq\n'%paths.avl_out
        ps.stdin.write(cmd.encode('utf-8'))
        ps.communicate()
        ps.kill()

        raw_data = ''
        fid = open(paths.avl_out,'rt')
        for line in fid.readlines():
            raw_data += line + ' '
        fid.close()
        raw_data.replace('\n',' ')
        idx = raw_data.find('CLtot =')
        CL = float(raw_data[idx+8:idx+18])
        idx = raw_data.find('CDtot =')
        CD = float(raw_data[idx+8:idx+18])
        idx = raw_data.find('Cmtot =')
        CM = float(raw_data[idx+8:idx+18])
        idx = raw_data.find('Xnp =')
        NP = float(raw_data[idx+8:idx+18])
        SM = (NP - self.cb.aircraft_cgx)/self.cb.wing_mac
        os.remove(paths.avl_out)
        return CL, CD, CM, SM

    def calc_avl_trim(self, velocity):
        q = 0.5*self.cb.air_density * velocity*velocity
        CL_req = self.cb.mass*self.cb.grav_accel / (self.cb.wing_area*q)
        ps=sp.Popen(paths.avl, stdin=sp.PIPE, stdout=sp.PIPE, 
                    stderr=sp.PIPE)
        cmd = 'load  %s\n'%paths.avl_inp  
        cmd += 'oper\n'
        cmd += 'm\nmn\n%.6f\n'%(velocity/self.cb.speed_of_sound)
        cmd += 'v\n%.6f\n'%self.cb.velocity
        cmd += 'd\n%.6f\n'%self.cb.air_density
        cmd += 'g\n%.6f\n'%self.cb.grav_accel
        cmd += 'm\n%.6f\n\n'%self.cb.mass
        cmd += 'a\nc\n%.6f\n'%CL_req
        cmd += 'd1\npm\n0\n'
        cmd += 'x\nst\n%s\n\nq\n'%paths.avl_out
        ps.stdin.write(cmd.encode('utf-8'))
        ps.communicate()
        ps.kill()

        raw_data = ''
        fid = open(paths.avl_out,'rt')
        for line in fid.readlines():
            raw_data += line + ' '
        fid.close()
        os.remove(paths.avl_out)
        raw_data.replace('\n',' ')
        idx = raw_data.find('Alpha =')
        aoa = float(raw_data[idx+8:idx+18])
        idx = raw_data.find('CLtot =')
        CL = float(raw_data[idx+8:idx+18])
        idx = raw_data.find('elevator        =')
        elevator = float(raw_data[idx+20:idx+26])
        return aoa, CL, elevator

    def run(self):
        self.set_flight_condition()
        self.calc_friction_drag()
        self.create_avl_input_file()

        fid = open(sys.argv[2]+'aero.out','wt')
        fid.write('aoa\tCL\tCD\tCM\tSM\n')
        aoa_list = np.arange(self.cb.aoa_min, self.cb.aoa_max, 
            self.cb.aoa_step)

        for aoa in aoa_list:
            cl, cd, cm, SM = self.calc_avl_aero(aoa)
            fid.write('%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n'%(aoa, cl, cd, 
                                                        cm, SM))
        fid.close()

        fid = open(sys.argv[2]+'trim.out','wt')
        fid.write('velocity\taoa\tCL\televator\n')
        velocity_list = np.arange(self.cb.velocity_min, 
            self.cb.velocity_max, self.cb.velocity_step)
        for velocity in velocity_list:
            aoa, CL, elevator = self.calc_avl_trim(velocity)
            fid.write('%.6f\t%.6f\t%.6f\t%.6f\n'%(velocity, aoa, CL,
                                                   elevator))
        fid.close()
        os.remove(paths.avl_inp)
        os.remove(paths.avl_fus_profile)

        data = np.genfromtxt(sys.argv[2]+'trim.out', skip_header=1)
        plt.figure()
        plt.grid()
        plt.plot(data[:,0], data[:,1], 'b-')
        plt.plot(data[:,0], data[:, 3], 'r-')
        plt.legend(['AoA', 'Trim elevator defl'])
        plt.xlabel('Velocity, m/s')
        plt.ylabel('Angle, deg')
        plt.savefig(sys.argv[2]+'trim.png')

if __name__=="__main__":
    ia = IntegratedAnalysis()
    ia.run()