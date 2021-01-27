#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Mikhail Zhigun"

import sys

assert sys.version_info[0] >= 3 and sys.version_info[1] >= 6, 'Python >= 3.6 is required'

from enum import Enum, IntEnum
from typing import NamedTuple, Union, Tuple
from abc import abstractmethod


class FPPrecision(IntEnum):
    SINGLE = 1
    DOUBLE = 2


class Real(NamedTuple):
    fp_precision: FPPrecision

    @property
    def size(self) -> int:
        if self.fp_precision == FPPrecision.SINGLE:
            return 4
        else:
            return 8


class Logical:
    @property
    def size(self) -> int:
        return 4


class Integer:
    @property
    def size(self) -> int:
        return 4


KB = 1024
MB = 1024 * 1024
GB = 1024 * 1024 * 1024


def size_str(size: int) -> str:
    if size < KB:
        return f'{size} bytes'
    elif size < MB:
        return f'{round(float(size) / KB, 2)} kb'
    elif size < GB:
        return f'{round(float(size) / MB, 2)} mb'
    else:
        return f'{round(float(size) / GB, 2)} gb'


def print_breakdown(val,  offset: int=0):
    print('\t' * offset + f'{val.name if hasattr(val, "name") else type(val).__name__} : {size_str(val.size)}')
    if isinstance(val, UserType):
        for member in val.members:
            print_breakdown(member, offset + 1)
    elif isinstance(val, Scalar):
        if isinstance(val.type, UserType):
            for member in val.type.members:
                print_breakdown(member, offset + 1)


class Scalar(NamedTuple):
    name: str
    type: Union[Real, Integer, Logical]

    @property
    def size(self) -> int:
        return self.type.size


class Dim(NamedTuple):
    name: str
    length: int


class Array(NamedTuple):
    name: str
    type: Union[Real, Integer, Logical]
    dims: Tuple[Dim]

    @property
    def size(self) -> int:
        sz = self.type.size
        for dim in self.dims:
            sz = sz * dim.length
        return sz


class UserType:
    @property
    def members(self): return self._members

    @property
    def size(self) -> int:
        sz = 0
        for member in self.members:
            sz = sz + member.size
        return sz

    def print_breakdown(self, offset=0):
        print_breakdown(self)

    def __init__(self):
        seld._members = []


class Single_level_type(UserType):
    def __init__(self, fp_precision: FPPrecision, ncol: int, nalbedobands: int, nemissbands: int, n_bands_sw: int):
        Cols = Dim('cols', ncol)
        Albedos = Dim('albedos', nalbedobands)
        Emissbands = Dim('emissbands', nemissbands)
        BandsSW = Dim('bands_sw', n_bands_sw)
        real = Real(fp_precision)
        int = Integer()
        self._members = (Array('cos_sza', real, (Cols,)),
                         Array('sw_albedo', real, (Cols, Albedos)),
                         Array('sw_albedo_direct', real, (Cols, Albedos)),
                         Array('lw_emissivity', real, (Cols, Emissbands)),
                         Array('lw_emission', real, (Cols, Emissbands)),
                         Array('solar_irradiance', real, (BandsSW,)),
                         Scalar('iseed', int),
                         Scalar('is_simple_surface', Logical()))


class Thermodynamics_type(UserType):
    def __init__(self, fp_precision: FPPrecision, ncol: int, nlev: int):
        Cols = Dim('cols', ncol)
        Levels = Dim('levels', nlev)
        Levels1 = Dim('levels1', nlev + 1)
        real = Real(fp_precision)
        self._members = (Array('pressure_hl', real, (Levels1,)),
                         Array('temperature_hl', real, (Levels1,)),
                         Array('h2o_sat_liq', real, (Cols, Levels)))


class Gas_type(UserType):
    def __init__(self, fp_precision: FPPrecision, ncol: int, nlev: int, NMaxGases: int):
        Cols = Dim('cols', ncol)
        Levels = Dim('levels', nlev)
        Gases = Dim('gases', NMaxGases)
        real = Real(fp_precision)
        integer = Integer()
        self._members = (Array('iunits', integer, (Gases,)),
                         Array('scale_factor', real, (Gases,)),
                         Array('mixing_ratio', real, (Cols, Levels, Gases)),
                         Array('is_present', Logical(), (Gases,)),
                         Array('is_well_mixed', Logical(), (Gases,)),
                         Scalar('ntype', integer),
                         Scalar('ncol', integer),
                         Scalar('nlev', integer),
                         Array('icode', integer, (Gases,)))


class Aerosol_type(UserType):
    def __init__(self, fp_precision: FPPrecision, ncol: int, nlev: int, n_aerosol_types: int, is_direct: bool,
                 n_bands_sw: int, n_bands_lw: int):
        Cols = Dim('cols', ncol)
        Levels = Dim('levels', nlev)
        SWBands = (Dim('sw_bands', n_bands_sw if is_direct else 0), Levels, Cols)
        LWBands = (Dim('lw_bands', n_bands_lw if is_direct else 0), Levels, Cols)
        AerosolTypes = Dim('aerosol_type', n_aerosol_types)
        real = Real(fp_precision)
        integer = Integer()
        self._members = (Array('mixing_ratio', real, (Cols, Levels, AerosolTypes)),
                         Array('od_sw', real, SWBands),
                         Array('ssa_sw', real, SWBands),
                         Array('g_sw', real, SWBands),
                         Array('od_lw', real, LWBands),
                         Array('ssa_lw', real, LWBands),
                         Array('g_lw', real, LWBands),
                         Scalar('istartlev', integer),
                         Scalar('iendlev', integer),
                         Scalar('is_direct', Logical()))


class Cloud_type(UserType):
    def __init__(self, fp_precision: FPPrecision, ncol: int, nlev: int):
        Cols = Dim('cols', ncol)
        Levels = Dim('levels', nlev)
        dim = (Cols, Levels)
        real = Real(fp_precision)
        self._members = (Array('q_liq', real, dim),
                         Array('q_ice', real, dim),
                         Array('re_liq', real, dim),
                         Array('re_ice', real, dim),
                         Array('fraction', real, dim),
                         Array('fractional_std', real, dim),
                         Array('inv_cloud_effective_size', real, dim),
                         Array('inv_inhom_effective_size', real, dim),
                         Array('overlap_param', real, dim))


class Flux_type(UserType):
    def __init__(self, fp_precision: FPPrecision, ncol: int, nlev: int, do_save_spectral_flux: bool, do_clear: bool,
                 do_lw: bool, n_spec_lw: int, do_lw_derivatives: bool, n_g_lw: int,
                 do_canopy_fluxes_lw: bool, n_canopy_bands_lw: int,
                 do_sw: bool, do_sw_direct: bool, n_spec_sw: int, do_surface_sw_spectral_flux: bool, n_bands_sw: int,
                 n_g_sw: int,
                 do_canopy_fluxes_sw: bool, n_canopy_bands_sw: int):
        Cols = Dim('cols', ncol)
        Levels = Dim('levels', nlev + 1)
        dimColLev = (Cols, Levels)
        real = Real(fp_precision)
        dimLWBands = Dim('lw_bands', n_spec_lw)
        dimSWBands = Dim('sw_bands', n_spec_sw)
        dimG_lw = Dim("g_lw", n_g_lw)
        dimG_sw = Dim("g_sw", n_g_sw)
        members = []
        if do_lw:
            members.append(Array('lw_up', real, dimColLev))
            members.append(Array('lw_dn', real, dimColLev))
            if do_clear:
                members.append(Array('lw_up_clear', real, dimColLev))
                members.append(Array('lw_dn_clear', real, dimColLev))
            if do_save_spectral_flux:
                members.append(Array('lw_up_band', real, (dimLWBands, Cols, Levels)))
                members.append(Array('lw_dn_band', real, (dimLWBands, Cols, Levels)))
                if do_clear:
                    members.append(Array('lw_up_clear_band', real, (dimLWBands, Cols, Levels)))
                    members.append(Array('lw_dn_clear_band', real, (dimLWBands, Cols, Levels)))
            if do_lw_derivatives:
                members.append(Array('lw_derivatives', real, (Cols, Levels)))
            members.append(Array('lw_dn_surf_g', real, (dimG_lw, Cols)))
            if do_clear:
                members.append(Array('lw_dn_surf_clear_g', real, (dimG_lw, Cols)))
            if do_canopy_fluxes_lw:
                members.append(Array('lw_dn_surf_canopy', real, (Dim("canopy_bands_lw", n_canopy_bands_lw), Cols)))
        if do_sw:
            members.append(Array('sw_up', real, dimColLev))
            members.append(Array('sw_dn', real, dimColLev))
            if do_sw_direct:
                members.append(Array('sw_dn_direct', real, dimColLev))
            if do_clear:
                members.append(Array('sw_up_clear', real, dimColLev))
                members.append(Array('sw_dn_clear', real, dimColLev))
                if do_sw_direct:
                    members.append(Array('sw_dn_direct_clear', real, dimColLev))
            if do_save_spectral_flux:
                members.append(Array('sw_up_band', real, (dimSWBands, Cols, Levels)))
                members.append(Array('sw_dn_band', real, (dimSWBands, Cols, Levels)))
                if do_sw_direct:
                    members.append(Array('sw_dn_direct_band', real, (dimSWBands, Cols, Levels)))
                if do_clear:
                    members.append(Array('sw_up_clear_band', real, (dimSWBands, Cols, Levels)))
                    members.append(Array('sw_dn_clear_band', real, (dimSWBands, Cols, Levels)))
                    if do_sw_direct:
                        members.append(Array('sw_dn_direct_clear_band', real, (dimSWBands, Cols, Levels)))
            if do_surface_sw_spectral_flux:
                dimSurfSWBands = Dim('sw_surf_bands', n_bands_sw)
                members.append(Array('sw_dn_surf_band', real, (dimSurfSWBands, Cols)))
                members.append(Array('sw_dn_direct_surf_band', real, (dimSurfSWBands, Cols)))
                if do_clear:
                    members.append(Array('sw_dn_surf_clear_band', real, (dimSurfSWBands, Cols)))
                    members.append(Array('sw_dn_direct_surf_clear_band', real, (dimSurfSWBands, Cols)))
            members.append(Array('sw_dn_diffuse_surf_g', real, (dimG_sw, Cols)))
            members.append(Array('sw_dn_direct_surf_g', real, (dimG_sw, Cols)))
            if do_clear:
                members.append(Array('sw_dn_diffuse_surf_clear_g', real, (dimG_sw, Cols)))
                members.append(Array('sw_dn_direct_surf_clear_g', real, (dimG_sw, Cols)))
            if do_canopy_fluxes_sw:
                members.append(
                    Array('sw_dn_diffuse_surf_canopy', real, (Dim("canopy_bands_sw", n_canopy_bands_sw), Cols)))
                members.append(
                    Array('sw_dn_direct_surf_canopy', real, (Dim("canopy_bands_sw", n_canopy_bands_sw), Cols)))
            members.append(Array('cloud_cover_lw', real, (Cols,)))
            members.append(Array('cloud_cover_sw', real, (Cols,)))
        self._members = tuple(members)


class RadiationIO(UserType):
    def __init__(self, fp_precision: FPPrecision, ncol: int, nlev: int, do_save_spectral_flux: bool, do_clear: bool,
                 do_lw: bool, n_spec_lw: int, do_lw_derivatives: bool, n_g_lw: int,
                 do_canopy_fluxes_lw: bool, n_canopy_bands_lw: int, n_bands_lw: int,
                 do_sw: bool, do_sw_direct: bool, n_spec_sw: int, do_surface_sw_spectral_flux: bool, n_bands_sw: int,
                 n_g_sw: int,
                 do_canopy_fluxes_sw: bool, n_canopy_bands_sw: int,
                 NMaxGases: int, n_aerosol_types: int, aerosol_is_direct: bool,
                 nalbedobands: int, nemissbands: int):
        self._members = (Scalar('ncol', Integer()),
                         Scalar('nlev', Integer()),
                         Scalar('istartcol', Integer()),
                         Scalar('iendcol', Integer()),
                         Scalar('single_level',
                                Single_level_type(fp_precision, ncol, nalbedobands, nemissbands, n_bands_sw)),
                         Scalar('thermodynamics', Thermodynamics_type(fp_precision, ncol, nlev)),
                         Scalar('gas', Gas_type(fp_precision, ncol, nlev, NMaxGases)),
                         Scalar('aerosol',
                                Aerosol_type(fp_precision, ncol, nlev, n_aerosol_types, aerosol_is_direct, n_bands_sw,
                                             n_bands_lw)),
                         Scalar('cloud', Cloud_type(fp_precision, ncol, nlev)),
                         Scalar('flux', Flux_type(fp_precision, ncol, nlev, do_save_spectral_flux, do_clear,
                                                  do_lw, n_spec_lw, do_lw_derivatives, n_g_lw,
                                                  do_canopy_fluxes_lw, n_canopy_bands_lw,
                                                  do_sw, do_sw_direct, n_spec_sw, do_surface_sw_spectral_flux,
                                                  n_bands_sw,
                                                  n_g_sw,
                                                  do_canopy_fluxes_sw, n_canopy_bands_sw)))


DEFAULT_CFG = {'fp_precision': FPPrecision.DOUBLE, 'ncol': 20000, 'nlev': 137, 'do_save_spectral_flux': False,
               'do_clear': True,
               'do_lw': True, 'n_spec_lw': 0, 'do_lw_derivatives': False, 'n_g_lw': 140, 'n_bands_lw': 16,
               'do_canopy_fluxes_lw': False, 'n_canopy_bands_lw': 1,
               'do_sw': True, 'do_sw_direct': True, 'n_spec_sw': 0, 'do_surface_sw_spectral_flux': False,
               'n_bands_sw': 14,
               'n_g_sw': 112, 'do_canopy_fluxes_sw': False, 'n_canopy_bands_sw': 1,
               'NMaxGases': 12, 'n_aerosol_types': 11, 'aerosol_is_direct': False,
               'nalbedobands': 2, 'nemissbands': 1}

if __name__ == '__main__':
    print_breakdown(RadiationIO(**DEFAULT_CFG))
    # print_breakdown(Scalar('bla', Integer()))
    # print_breakdown(Array('arr', Real(FPPrecision.DOUBLE), (Dim('Cols', 10), Dim('Rows', 10), Dim('Lvls', 20))))
    # print(RadiationIO(**DEFAULT_CFG).size)
    # print(Real(FPPrecision.SINGLE).size)
    # print(Real(FPPrecision.DOUBLE).size)
    # print(Integer().size)
    # print(Logical().size)
    # print(Scalar('bla', Integer()).size)
    # print(Scalar('bla', Logical()).size)
    # print(Array('arr', Real(FPPrecision.DOUBLE), (Dim('Cols', 10), Dim('Rows', 10), Dim('Lvls', 20))).size)
    # slt = Single_level_type(FPPrecision.DOUBLE, 32, 2, 1, 140)
    # print(slt.size)
