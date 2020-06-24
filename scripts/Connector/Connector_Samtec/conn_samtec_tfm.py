#!/usr/bin/env python

'''
kicad-footprint-generator is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

kicad-footprint-generator is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with kicad-footprint-generator. If not, see < http://www.gnu.org/licenses/ >.
'''

import sys
import os
#sys.path.append(os.path.join(sys.path[0],"..","..","kicad_mod")) # load kicad_mod path

# export PYTHONPATH="${PYTHONPATH}<path to kicad-footprint-generator directory>"
sys.path.append(os.path.join(sys.path[0], "..", "..", ".."))  # load parent path of KicadModTree
from math import sqrt
import argparse
import yaml
from helpers import *
from KicadModTree import *

sys.path.append(os.path.join(sys.path[0], "..", "..", "tools"))  # load parent path of tools
from footprint_text_fields import addTextFields

series = 'TFM'
series_long = 'High-Reliability Tiger Eye Terminal Strips, .050" Pitch'
manufacturer = 'Samtec'
orientation = 'H'
number_of_rows = 2
pitch = 1.27

datasheet_tht = 'http://suddendocs.samtec.com/prints/tfm-1xx-xx-xxx-d-xxx-mkt.pdf'
footprint_tht = 'http://suddendocs.samtec.com/prints/tfm(x)-1xx-x1-xx-d-xx-footprint.pdf'

# BE = Bottem Enty
# TE = Top Entry
# PE = Pass Through Entry

drill = 0.635

pad_to_pad_clearance = 1.2
max_annular_ring = 0.3
min_annular_ring = 0.15

variant_params = {
    'tht-a': {
        'pins_per_row_range': [15, 20],
        'pad_drill': drill,
        'datasheets': [datasheet_tht, footprint_tht],
        'part_code': "TFM-1{n:02}-01-F-D-A",
        'peg_drill': 1.42,
        },
}

def generate_one_footprint(pins_per_row, variant, configuration):
    is_smd = variant_params[variant].get('smd', False)

    mpn = variant_params[variant]['part_code'].format(n=pins_per_row)
    alt_mpn = [code.format(n=pins_per_row) for code in variant_params[variant].get('alternative_codes', [])]

    # handle arguments
    orientation_str = configuration['orientation_options'][orientation]
    footprint_name = configuration['fp_name_format_string'].format(
        man=manufacturer,
        series='',
        mpn=mpn,
        num_rows=number_of_rows,
        pins_per_row=pins_per_row,
        pitch=pitch,
        mounting_pad = "",
        orientation=orientation_str)
    footprint_name = footprint_name.replace('__','_')

    kicad_mod = Footprint(footprint_name)
    kicad_mod.setDescription("{manufacturer} {series}, {mpn}{alt_mpn}, {pins_per_row} Pins per row ({datasheet}), generated with kicad-footprint-generator".format(
        manufacturer = manufacturer,
        series = series_long,
        mpn = mpn,
        alt_mpn = ' (compatible alternatives: {})'.format(', '.join(alt_mpn)) if len(alt_mpn) > 0 else "",
        pins_per_row = pins_per_row,
        datasheet = ', '.join(variant_params[variant]['datasheets'])))

    kicad_mod.setTags(configuration['keyword_fp_string'].format(series=series,
        orientation=orientation_str, man=manufacturer,
        entry=configuration['entry_direction'][orientation]))

    if is_smd:
        kicad_mod.setAttribute('smd')

    ########################## Dimensions ##############################

    # Solder pads pitch
    pad_pitch_x = pitch
    pad_pitch_y = variant_params[variant].get('pad_pitch_y', -pitch)

    # Size of the connector (without pins)
    size_x = pins_per_row * pitch + 3.18
    size_y = 5.72

    # Offset everything
    offset_x = 0
    offset_y = (pad_pitch_y - pitch)/2

    # Position of connector pins. This is either used for npth holes at pin position or for tht pads
    pin_row1_y = offset_y + 0
    pin_row2_y = pin_row1_y + pitch
    pin1_x = offset_x + 0

    # By definition pad1 needs to be at 0,0 for THT
    pad1_x = 0
    pad_row1_y = 0
    pad_row2_y = pad_row1_y + pad_pitch_y
    pad_layer = Pad.LAYERS_THT
    pad_drill = variant_params[variant].get('pad_drill', None)
    pad_type = Pad.TYPE_THT

    pad_size = [pitch - pad_to_pad_clearance, pad_drill + 2*max_annular_ring]
    if pad_size[0] - pad_drill < 2*min_annular_ring:
        pad_size[0] = pad_drill + 2*min_annular_ring
    if pad_size[0] - pad_drill > 2*max_annular_ring:
        pad_size[0] = pad_drill + 2*max_annular_ring

    if pad_size[1] - pad_drill < 2*min_annular_ring:
        pad_size[1] = pad_drill + 2*min_annular_ring
    if pad_size[1] - pad_drill > 2*max_annular_ring:
        pad_size[1] = pad_drill + 2*max_annular_ring

    pad_shape=Pad.SHAPE_CIRCLE

    peg_distance = pins_per_row * pitch + 1.91
    peg1_x = offset_x + pins_per_row/2 * pitch - peg_distance/2 - pitch/2
    peg2_x = offset_x + pins_per_row/2 * pitch + peg_distance/2 - pitch/2
    peg_y = offset_y + pitch/2

    body_edge = {
        'left': offset_x + pins_per_row/2 * pitch - size_x/2 - pitch/2,
        'right': offset_x + pins_per_row/2 * pitch + size_x/2 - pitch/2,
        'top': offset_y - size_y/2 + pitch/2,
        'bottom': offset_y + size_y/2 + pitch/2,
    }

    bounding_box={
        'left': body_edge['left'],
        'right': body_edge['right'],
        'top': min(pad_row1_y - pad_size[1]/2, body_edge['top']),
        'bottom': max(pad_row2_y + pad_size[1]/2, body_edge['bottom']),
        }

    # ############################# Pads ##################################

    # Pegs
    peg_drill = variant_params[variant].get('peg_drill', None)
    if peg_drill != None:
        kicad_mod.append(Pad(at=[peg1_x, peg_y], number="",
            type=Pad.TYPE_NPTH, shape=Pad.SHAPE_CIRCLE, size=peg_drill,
            drill=peg_drill, layers=Pad.LAYERS_NPTH))
        if (peg1_x != peg2_x):
            kicad_mod.append(Pad(at=[peg2_x, peg_y], number="",
                type=Pad.TYPE_NPTH, shape=Pad.SHAPE_CIRCLE, size=peg_drill,
                drill=peg_drill, layers=Pad.LAYERS_NPTH))

    # NPTH Holes for Connector Pins
    npth_drill = variant_params[variant].get('npth_drill', None)
    if npth_drill:
        for x in range(0, pins_per_row):
            for y in range(0, number_of_rows):
                kicad_mod.append(Pad(at=[pin1_x+x*pitch, pin_row1_y+y*pitch], number="",
                    type=Pad.TYPE_NPTH, shape=Pad.SHAPE_CIRCLE, size=npth_drill,
                    drill=npth_drill, layers=Pad.LAYERS_NPTH))

    # Pads
    optional_pad_params = {}
    if not is_smd:
        if configuration['kicad4_compatible']:
            optional_pad_params['tht_pad1_shape'] = Pad.SHAPE_RECT
        else:
            optional_pad_params['tht_pad1_shape'] = Pad.SHAPE_CIRCLE

    kicad_mod.append(PadArray(start=[pad1_x, pad_row1_y], initial=1,
        pincount=pins_per_row, increment=2,  x_spacing=pitch, size=pad_size,
        type=pad_type, shape=pad_shape, layers=pad_layer, drill=pad_drill,
        **optional_pad_params))
    kicad_mod.append(PadArray(start=[pad1_x, pad_row2_y], initial=2,
        pincount=pins_per_row, increment=2, x_spacing=pitch, size=pad_size,
        type=pad_type, shape=pad_shape, layers=pad_layer, drill=pad_drill,
        **optional_pad_params))

    # ######################## Fabrication Layer ###########################

    poly_f = [
        {'x': body_edge['left'], 'y': body_edge['top']},
        {'x': body_edge['right'], 'y': body_edge['top']},
        {'x': body_edge['right'], 'y': body_edge['bottom']},
        {'x': body_edge['left'], 'y': body_edge['bottom']},
        {'x': body_edge['left'], 'y': body_edge['top']},
    ]
    kicad_mod.append(PolygoneLine(polygone=poly_f,
        width=configuration['fab_line_width'], layer="F.Fab"))

    # Pin 1 marking
    p1m_sl = 1
    p1m_poly = [
        {'x': pad1_x - p1m_sl/2, 'y': body_edge['top']},
        {'x': pad1_x, 'y': body_edge['top'] + p1m_sl/sqrt(2)},
        {'x': pad1_x + p1m_sl/2, 'y': body_edge['top']}
    ]
    kicad_mod.append(PolygoneLine(polygone=p1m_poly,
        width=configuration['fab_line_width'], layer="F.Fab"))

    ############################ SilkS ##################################

    s_pad_offset = configuration['silk_pad_clearance'] + configuration['silk_line_width']/2
    s_offset = configuration['silk_fab_offset']

    # Check if we can draw a solid box for the connector or we need to interrupt it for pads
    if ((pad_row1_y - pad_size[1]/2 - s_pad_offset) < body_edge['top']) and ((pad_row1_y + pad_size[1]/2 + s_pad_offset) > body_edge['top']):
        s_xp1_left = pad1_x - pad_size[0]/2 - s_pad_offset
        s_xpn_right = pad1_x + size_x - pitch + pad_size[0]/2 + s_pad_offset

        poly_s_l = [
            {'x': s_xp1_left, 'y': pad_row1_y - pad_size[1]/2},
            {'x': s_xp1_left, 'y': body_edge['top'] - s_offset},
            {'x': body_edge['left'] - s_offset, 'y': body_edge['top'] - s_offset},
            {'x': body_edge['left'] - s_offset, 'y': body_edge['bottom'] + s_offset},
            {'x': s_xp1_left, 'y': body_edge['bottom'] + s_offset},
        ]
        kicad_mod.append(PolygoneLine(polygone=poly_s_l,
            width=configuration['silk_line_width'], layer="F.SilkS"))

        poly_s_r = [
            {'x': s_xpn_right, 'y': body_edge['top'] - s_offset},
            {'x': body_edge['right'] + s_offset, 'y': body_edge['top'] - s_offset},
            {'x': body_edge['right'] + s_offset, 'y': body_edge['bottom'] + s_offset},
            {'x': s_xpn_right, 'y': body_edge['bottom'] + s_offset},
        ]
        kicad_mod.append(PolygoneLine(polygone=poly_s_r,
            width=configuration['silk_line_width'], layer="F.SilkS"))

        for i in range(0, pins_per_row-1):
            s_xpi_left = pad1_x + i * pitch + pad_size[0]/2 + s_pad_offset
            s_xpi_right = pad1_x + (i+1) * pitch - pad_size[0]/2 - s_pad_offset

            poly_s_t = [
                {'x': s_xpi_left, 'y': body_edge['top'] - s_offset},
                {'x': s_xpi_right, 'y': body_edge['top'] - s_offset},
            ]
            kicad_mod.append(PolygoneLine(polygone=poly_s_t,
                width=configuration['silk_line_width'], layer="F.SilkS"))

            poly_s_b = [
                {'x': s_xpi_left, 'y': body_edge['bottom'] + s_offset},
                {'x': s_xpi_right, 'y': body_edge['bottom'] + s_offset},
            ]
            kicad_mod.append(PolygoneLine(polygone=poly_s_b,
                width=configuration['silk_line_width'], layer="F.SilkS"))

    else:
        poly_s = [
            {'x': body_edge['left'] - s_offset, 'y': body_edge['top'] - s_offset},
            {'x': body_edge['right'] + s_offset, 'y': body_edge['top'] - s_offset},
            {'x': body_edge['right'] + s_offset, 'y': body_edge['bottom'] + s_offset},
            {'x': body_edge['left'] - s_offset, 'y': body_edge['bottom'] + s_offset},
            {'x': body_edge['left'] - s_offset, 'y': body_edge['top'] - s_offset},
        ]
        kicad_mod.append(PolygoneLine(polygone=poly_s,
            width=configuration['silk_line_width'], layer="F.SilkS"))

    # Draw inner silkscreen showing notch in connector
    s_gap_height = 1.91
    s_notch_width = 2.79
    s_notch_height = 0.51
    s_inside_height = 3.43
    s_inside_width = pins_per_row * pitch + 0.64
    center_x = offset_x + pins_per_row/2 * pitch - pitch/2
    center_y = pad_pitch_y/2
    s_inside_left_x = center_x - s_inside_width/2
    s_inside_right_x = center_x + s_inside_width/2
    s_inside_top_y = center_y - s_inside_height/2
    s_inside_bottom_y = center_y + s_inside_height/2
    poly_s_inside_t = [
        {'x': s_inside_left_x, 'y': s_inside_top_y + s_gap_height/2},
        {'x': s_inside_left_x, 'y': s_inside_top_y},
        {'x': s_inside_right_x - s_notch_width, 'y': s_inside_top_y},
        {'x': s_inside_right_x - s_notch_width, 'y': s_inside_top_y - s_notch_height},
        {'x': s_inside_right_x, 'y': s_inside_top_y - s_notch_height},
        {'x': s_inside_right_x, 'y': s_inside_top_y + s_gap_height/2}
    ]
    kicad_mod.append(PolygoneLine(polygone=poly_s_inside_t,
        width=configuration['silk_line_width'], layer="F.SilkS"))

    poly_s_inside_b = [
        {'x': s_inside_left_x, 'y': s_inside_bottom_y - s_gap_height/2},
        {'x': s_inside_left_x, 'y': s_inside_bottom_y},
        {'x': s_inside_right_x, 'y': s_inside_bottom_y},
        {'x': s_inside_right_x, 'y': s_inside_bottom_y - s_gap_height/2}
    ]
    kicad_mod.append(PolygoneLine(polygone=poly_s_inside_b,
        width=configuration['silk_line_width'], layer="F.SilkS"))


    # ############################ CrtYd ##################################

    cy_offset = configuration['courtyard_offset']['connector']
    cy_grid = configuration['courtyard_grid']

    cy_top = roundToBase(bounding_box['top'] - cy_offset, cy_grid)
    cy_bottom = roundToBase(bounding_box['bottom'] + cy_offset, cy_grid)
    cy_left = roundToBase(bounding_box['left'] - cy_offset, cy_grid)
    cy_right = roundToBase(bounding_box['right'] + cy_offset, cy_grid)

    poly_cy = [
        {'x': cy_left, 'y': cy_top},
        {'x': cy_right, 'y': cy_top},
        {'x': cy_right, 'y': cy_bottom},
        {'x': cy_left, 'y': cy_bottom},
        {'x': cy_left, 'y': cy_top},
    ]
    kicad_mod.append(PolygoneLine(polygone=poly_cy,
        layer='F.CrtYd', width=configuration['courtyard_line_width']))

    ######################### Text Fields ###############################

    addTextFields(kicad_mod=kicad_mod, configuration=configuration, body_edges=body_edge,
        courtyard={'top':cy_top, 'bottom':cy_bottom}, fp_name=footprint_name, text_y_inside_position='bottom')

    ##################### Output and 3d model ############################

    model3d_path_prefix = configuration.get('3d_model_prefix','${KISYS3DMOD}/')
    lib_name_suffix = '_SMD' if is_smd else '_THT'

    lib_name = configuration['lib_name_format_string_full'].format(series=series, man=manufacturer, suffix=lib_name_suffix)

    model_name = '{model3d_path_prefix:s}{lib_name:s}.3dshapes/{fp_name:s}.wrl'.format(
        model3d_path_prefix=model3d_path_prefix, lib_name=lib_name, fp_name=footprint_name)
    kicad_mod.append(Model(filename=model_name))

    output_dir = '{lib_name:s}.pretty/'.format(lib_name=lib_name)
    if not os.path.isdir(output_dir): #returns false if path does not yet exist!! (Does not check path validity)
        os.makedirs(output_dir)
    filename =  '{outdir:s}{fp_name:s}.kicad_mod'.format(outdir=output_dir, fp_name=footprint_name)

    file_handler = KicadFileHandler(kicad_mod)
    file_handler.writeFile(filename)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='use confing .yaml files to create footprints.')
    parser.add_argument('--global_config', type=str, nargs='?', help='the config file defining how the footprint will look like. (KLC)', default='../../tools/global_config_files/config_KLCv3.0.yaml')
    parser.add_argument('--series_config', type=str, nargs='?', help='the config file defining series parameters.', default='../conn_config_KLCv3.yaml')
    parser.add_argument('--kicad4_compatible', action='store_true', help='Create footprints kicad 4 compatible')
    args = parser.parse_args()

    with open(args.global_config, 'r') as config_stream:
        try:
            configuration = yaml.safe_load(config_stream)
        except yaml.YAMLError as exc:
            print(exc)

    with open(args.series_config, 'r') as config_stream:
        try:
            configuration.update(yaml.safe_load(config_stream))
        except yaml.YAMLError as exc:
            print(exc)

    configuration['kicad4_compatible'] = args.kicad4_compatible

    for variant in variant_params:
        for pins_per_row in variant_params[variant]['pins_per_row_range']:
            generate_one_footprint(pins_per_row, variant, configuration)
