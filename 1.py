#!/usr/bin/env python3
#coding utf-8

import sys 
import math

import pprint
pprint = pprint.PrettyPrinter(indent=4).pprint

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

import lat
import lat.read_atoms
import lat.define_phase

SCALE = 20
RADIUS = 5
WIDTH = 5
BIG_INT = 1000000
SYSTEM_NAME = '10x20'
PHASE = 'modifier'
SCHEME_NAME = 'phase'
DRAW_HYDROGENS = True
CLIPPING_X = None
CLIPPING_Y = None
CLIPPING_Z = None

class Atom():
    def __init__(self, atomNumber,
                       atomGroup,
                       atomType,
                       atomCharge,
                       x, y, z):
        self.values = {}
        self.values['atomNumber'] = atomNumber
        self.values['atomGroup'] = atomGroup
        self.values['atomType'] = atomType
        self.values['atomCharge'] = atomCharge
        self.values['X'] = x
        self.values['Y'] = y
        self.values['Z'] = z

class Bounds():
    def __init__(self, xlo, xhi,
                       ylo, yhi,
                       zlo, zhi):
        self.values = {}
        self.values['Xlo'] = xlo
        self.values['Xhi'] = xhi
        self.values['Ylo'] = ylo
        self.values['Yhi'] = yhi
        self.values['Zlo'] = zlo
        self.values['Zhi'] = zhi
        self.values['lX'] = xhi - xlo
        self.values['lY'] = yhi - ylo
        self.values['lZ'] = zhi - zlo

class Bond():
    def __init__(self, bondNumber,
                       bondType,
                       bondAtomOneNumber,
                       bondAtomTwoNumber):
        self.values = {}
        self.values['bondNumber'] = bondNumber
        self.values['bondType'] = bondType
        self.values['bondAtomOneNumber'] = bondAtomOneNumber
        self.values['bondAtomTwoNumber'] = bondAtomTwoNumber
        self.values['bondBegin'] = {}
        self.values['bondEnd'] = {}

    def calculateEnds(self, atomsList):
        l = len(atomsList)
        atomOne = self.values['bondAtomOneNumber']
        atomTwo = self.values['bondAtomTwoNumber']
        self.values['bondBegin']['X'] = atomsList[atomOne].values['X']
        self.values['bondBegin']['Y'] = atomsList[atomOne].values['Y']
        self.values['bondBegin']['Z'] = atomsList[atomOne].values['Z']
        self.values['bondEnd']['X'] = atomsList[atomTwo].values['X']
        self.values['bondEnd']['Y'] = atomsList[atomTwo].values['Y']
        self.values['bondEnd']['Z'] = atomsList[atomTwo].values['Z']

class AtomSystem():
    def __init__(self, atomsList, bondsList, bounds):
        self.atomsList = atomsList
        self.bondsList = bondsList
        self.bounds = bounds

    def computeRanges(self):
        minX = minY = minZ = BIG_INT
        maxX = maxY = maxZ = -BIG_INT
        for atom in self.atomsList.values():
            if atom.values['X'] < minX:
                minX = atom.values['X']
            if atom.values['X'] > maxX:
                maxX = atom.values['X']
            if atom.values['Y'] < minY:
                minY = atom.values['Y']
            if atom.values['Y'] > maxY:
                maxY = atom.values['Y']
            if atom.values['Z'] < minZ:
                minZ = atom.values['Z']
            if atom.values['Z'] > maxZ:
                maxZ = atom.values['Z']
        self.ranges = {}
        self.ranges['minX'] = minX
        self.ranges['maxX'] = maxX
        self.ranges['minY'] = minY
        self.ranges['maxY'] = maxY
        self.ranges['minZ'] = minZ
        self.ranges['maxZ'] = maxZ

class MainWidget(QWidget):
    """Рисует атомы и связи"""
    def __init__(self, atomSystem, bounds,
                       coord1='X', coord2='Y'):
        super().__init__()
        self.atomSystem = atomSystem
        self.radius = RADIUS
        self.scale = SCALE
        self.width = WIDTH
        self.coord1 = coord1
        self.coord2 = coord2
        x = (self.scale * (atomSystem.ranges['max' + self.coord1] - 
                           atomSystem.ranges['min' + self.coord1]) +
            2 * self.radius + 
             self.width)
        y = (self.scale * (atomSystem.ranges['max' + self.coord2] -
                           atomSystem.ranges['min' + self.coord2]) +
             2 * self.radius + 
             self.width)
        self.resize(x, y)
        self.takeScreenshot(coord1 + coord2)

    def paintEvent(self, e):
        qp = QPainter()
        qp.begin(self)
        brush = QBrush(QColor(0, 0, 0), Qt.SolidPattern)
        pen = QPen(brush, WIDTH)
        qp.setPen(pen)
        dxyz = self.atomSystem.ranges
        for bond in self.atomSystem.bondsList.values():
            if SCHEME_NAME != 'full':
                continue
            start = bond.values['bondBegin']
            end = bond.values['bondEnd']
            length = ((start['X'] - end['X'])**2 + 
                      (start['Y'] - end['Y'])**2 +
                      (start['Z'] - end['Z'])**2)
            if length < 5:
                qp.drawLine(self.scale * (bond.values['bondBegin'][self.coord1] -
                                          dxyz['min' + self.coord1]) +
                            self.radius + 
                            self.width / 2,
                            self.scale * (bond.values['bondBegin'][self.coord2] -
                                          dxyz['min' + self.coord2]) +
                            self.radius + 
                            self.width / 2,
                            self.scale * (bond.values['bondEnd'][self.coord1] -
                                          dxyz['min' + self.coord1]) +
                            self.radius + 
                            self.width / 2,
                            self.scale * (bond.values['bondEnd'][self.coord2] -
                                          dxyz['min' + self.coord2]) +
                            self.radius + 
                            self.width / 2)
        l = len(self.atomSystem.atomsList.values())
        for atom in self.atomSystem.atomsList.values():
            if (CLIPPING_X is not None and
                (CLIPPING_X - self.radius > atom.values['X'] or
                 atom.values['X'] > CLIPPING_X + self.radius)):
                continue
            if (CLIPPING_Y is not None and
                (CLIPPING_Y - self.radius > atom.values['Y'] or
                 atom.values['Y'] > CLIPPING_Y + self.radius)):
                continue
            if (CLIPPING_Z is not None and
                (CLIPPING_Z - self.radius > atom.values['Z'] or
                 atom.values['Z'] > CLIPPING_Z + self.radius)):
                continue
            if SYSTEM_NAME == '10x20' and (SCHEME_NAME == 'full'
                                           or SCHEME_NAME == 'atomic'):
                if atom.values['atomType'] in [4, 5, 6, 7, 15]:
                    brush = QBrush(QColor(255, 0, 0), Qt.SolidPattern)
                elif atom.values['atomType'] in [9, 16, 17]:
                    brush = QBrush(QColor(0, 255, 0), Qt.SolidPattern)
                elif atom.values['atomType'] in [10, 11, 14]:
                    brush = QBrush(QColor(105, 105, 105), Qt.SolidPattern)
                elif atom.values['atomType'] in [1]:
                    brush = QBrush(QColor(255, 0, 255), Qt.SolidPattern)
                elif atom.values['atomType'] in [2]:
                    brush = QBrush(QColor(128, 0, 0), Qt.SolidPattern)
                elif atom.values['atomType'] in [3]:
                    brush = QBrush(QColor(205, 92, 92), Qt.SolidPattern)
                elif atom.values['atomType'] in [8, 12, 13]:
                    if DRAW_HYDROGENS == False:
                        continue
                    brush = QBrush(QColor(255, 255, 255), Qt.SolidPattern)
                else: # неописанные атомы заметного розового цвета
                    print(atom.values['atomType'])
                    brush = QBrush(QColor(255, 0, 255), Qt.SolidPattern)
            elif SYSTEM_NAME == '10x20' and SCHEME_NAME == 'phase':
                if lat.define_phase.define_phase(l,
                                                 atom.values['atomNumber']) == 1:
                    brush = QBrush(QColor(255, 0, 0), Qt.SolidPattern)
                elif lat.define_phase.define_phase(l,
                                                   atom.values['atomNumber']) == 2:
                    brush = QBrush(QColor(0, 255, 0), Qt.SolidPattern)
                elif lat.define_phase.define_phase(l,
                                                   atom.values['atomNumber']) == 3:
                    brush = QBrush(QColor(0, 0, 255), Qt.SolidPattern)
            qp.setBrush(brush)
            pen = QPen(brush, 5)
            qp.setPen(pen)
            qp.drawEllipse(QPoint(self.scale * (atom.values[self.coord1] -
                                                dxyz['min' + self.coord1]) +
                                  self.radius + 
                                  self.width / 2,
                                  self.scale * (atom.values[self.coord2] -
                                                dxyz['min' + self.coord2]) +
                                  self.radius + 
                                  self.width / 2),
                           self.radius, self.radius)
        qp.end()

    def takeScreenshot(self, fname='XY.png', format='png'):
        p = self.grab();
        fname2 = ('pics/' + 
                  self.coord1 +
                  self.coord2 +
                  '.png')
        p.save(fname2, format, -1)


def main():
    app = QApplication(sys.argv)

    fname = '/home/anton/DataExamples/10x20.data'
    if SYSTEM_NAME == '10x20':
        moleculeNumber = 90

    i = 0
    atomsList = {}
    bondsList = {}
    print('Parsing ', fname);
    [atoms, bounds, bonds, angles] = lat.read_atoms.read_atoms(fname)
    for atom in atoms:
        i += 1
        print('Preparing atoms: ', i, ' / ', len(bonds))
        at = Atom(atom[0], atom[1], atom[2],
                  atom[3], atom[4], atom[5], atom[6])
        atomsList[at.values['atomNumber']] = at
    bounds = Bounds(bounds[0], bounds[1], bounds[2],
                    bounds[3], bounds[4], bounds[5])
    i = 0
    for bond in bonds:
        i += 1
        print('Preparing bonds: ', i, ' / ', len(bonds))
        bo = Bond(bond[0], bond[1], bond[2], bond[3])
        bo.calculateEnds(atomsList)
        bondsList[bo.values['bondNumber']] = bo
    i = 0

    atomSystem = AtomSystem(atomsList, bondsList, bounds)
    atomSystem.computeRanges()

    w = MainWidget(atomSystem, bounds, 'X', 'Y')
    w = MainWidget(atomSystem, bounds, 'X', 'Z')
    w = MainWidget(atomSystem, bounds, 'Y', 'Z')

    app.quit()

main()
