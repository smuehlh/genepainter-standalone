module Template
    extend self

    def generate_color_exons_script(intron_pos, chain, grey_pos)
        script = first_part_general(intron_pos, chain, grey_pos)
        script << specific_for_exons
        script << last_part_general
        script << last_specific_for_exons
        return script
    end

    def generate_color_splicesites_script(intron_pos, chain, grey_pos)
        script = first_part_general(intron_pos, chain, grey_pos)
        script << specific_for_splicesites
        script << last_part_general
        script << last_specific_for_splicesites
        return script
    end

    def first_part_general(intron_pos, chain, grey_pos)
        str = <<-EOS
# Generated by GenePainter.
# Copyright (c) 2014 Stefanie Muehlhausen, Martin Kollmar

from pymol import cmd

intron_pos = {#{intron_pos}}
chain = "#{chain}"
no_data_pos = [#{grey_pos}]
color_no_info = "grey70"
        EOS
    end

    def specific_for_exons
        str = <<-EOS
def color_exons():
    """
    usage: color_exons

    Exons are colored.
    First exon: green.

    Residues w/o information about the gene structure: grey
    """

    cmd.hide("all")
    cmd.color(color_no_info, "all")

    colors = ['lime', 'olive', 'deepteal', 'marine', 'hotpink', 'grey50', 'yelloworange', 'violetpurple']
    ncolors = len(colors)

    pos_start = 1
    i_col = 0
    for i_pos in sorted(intron_pos):
        if i_pos <= pos_start:
            pos_start = i_pos
            continue
        sel = "resi %d-%d and chain %s" % (pos_start, i_pos, chain)
        use_color = colors[i_col]

        print sel, "->", use_color
        cmd.color(use_color,sel)
        pos_start = i_pos+1
        i_col = i_col+1
        if(i_col == ncolors):
            i_col = 0

    for i_pos in no_data_pos:
        sel = "resi %d and chain %s" % (i_pos, chain)
        print sel, "-> grey"
        cmd.color(color_no_info, sel)

        EOS
    end

    def specific_for_splicesites
        str = <<-EOS
def color_splicesite():
    """
    usage: color_splicesite 

    Intron phase:       Default color:
    phase 0             green - dark grenn
    phase 1             blue - dark blue
    phase 2             red - dark red

    Last residue of exon: ligth color
    First residue of next exon: dark color

    Residues w/o information about the gene structure: grey
    """

    colors = {
        0: "splitpea",
        1: "lightteal",
        2: "darksalmon",
        10: "forest",
        11: "deepblue",
        12: "raspberry",
        '?': "lightblue",
        '??': "slate"
    }

    cmd.hide("all")
    cmd.color(color_no_info, "all")

    cmd.set("label_font_id", "7")
    cmd.set("label_size", "20")
    cmd.set("label_outline_color", "black")
    cmd.set("label_color", "white")
    cmd.set("label_shadow_mode", "off")

    for i_pos in sorted(intron_pos):
        sel = "resi %d and chain %s" % (i_pos, chain)
        print sel, "->", colors[intron_pos[i_pos]]
        cmd.color(colors[intron_pos[i_pos]],sel)
        if (intron_pos[i_pos] < 10 or intron_pos[i_pos] == '?'):
            cmd.label("%s and name c" % sel,"'  %s'"% (intron_pos[i_pos]))

    for i_pos in no_data_pos:
        sel = "resi %d and chain %s" % (i_pos, chain)
        print sel, "->", color_no_info
        cmd.color(color_no_info,sel)

        EOS
    end

    def last_part_general
        str = <<-EOS
    cmd.show("cartoon", "all")
    cmd.set("cartoon_side_chain_helper", "on")
    cmd.set("cartoon_fancy_helices", "on")

    # chains of no interest
    sel = " not chain %s" % chain
    cmd.cartoon("loop", sel)

    cmd.set("specular", "off")
    cmd.set("antialias", "1")
    cmd.set("ray_shadows", "off")
    cmd.set("ray_trace_mode", "1")
        EOS
    end

    def last_specific_for_splicesites
        str = <<-EOS
cmd.extend("color_splicesites",color_splicesite)
        EOS
    end

    def last_specific_for_exons
        str = <<-EOS
cmd.extend("color_exons",color_exons)
        EOS
    end

    def successor_phase_for_nondescriptive_intron
        # important to have "untouched" quotes
        return "\"??\""
    end
end