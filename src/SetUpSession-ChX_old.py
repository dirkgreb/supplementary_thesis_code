from chimerax.core.commands import run


def SetUpSession(session):
    run(session,'set bgColor white; graphics silhouettes true; lighting soft depthCue false')
    run(session,"alias thin_licorice car style protein modeh default arrows f xsect oval width 0.5 thick 0.5")
    run(session,"alias licorice car style protein modeh default arrows f xsect oval width 1 thick 1")
    run(session,"alias thick_licorice car style protein modeh default arrows f xsect oval width 2 thick 2")
    run(session,"alias med_licorice car style protein modeh default arrows f xsect oval width 1.5 thick 1.5")
    run(session,"alias thin_licorice_chain car style $1 modeh default arrows f xsect oval width 0.75 thick 0.75")
    run(session,"alias save_png save $1.png supersample 4 transparentBackground true")
    run(session,"alias spin_pic save_png $1-1; turn y 90;  save_png $1-2 ; turn y 90;  save_png $1-3; turn y 90; save_png $1-4; turn y 90")
    run(session,"alias spin_pic_six save_png $1-1; turn y 60;  save_png $1-2 ; turn y 60;  save_png $1-3; turn y 60; save_png $1-4; turn y 60 ; save_png $1-5; turn y 60;  save_png $1-6 ; turn y 60")

    run(session,"alias match_af2s mmaker #2-5/$1 to #1/$1")
    
    run(session,"alias show_all_cartoon hide atoms ; hide surface ; show cartoon")
    
    run(session,"alias show_this_model hide cartoon ; hide atoms ; hide surface ; show $1 cartoon")
    
    run(session,"alias af2_color col bfactor $1 palette alphafold")
    
    run(session,"alias col_conservation color byattr seq_conservation $1 & protein palette cyanmaroon novalue black ; col $1 fromatoms target s")
    
    run(session,"alias col_conservation_bydomain color byattr seq_conservation $1/$2 & protein palette cyanmaroon novalue black ; col $1/$2 fromatoms target s")

    run(session,"alias show_conserved_atoms show /$1::seq_conservation>1.8 atoms ; show /$1 surface ; col fromatoms ; transparency 70 ")
    run(session,"alias cylinders preset cartoons/nucleotides cylinders/stubs")
    run(session,"alias ribbons preset cartoons/nucleotides ribbons/slabs")
    run(session,'alias mlp_bydomain mlp $1/$2') 
    run(session,'alias coul_bydomain coulombic $1/$2')
    
    run(session,'alias col_t7_doms_2 col $1/A #75c376 target acs; col $1/B #a7cee2 target acs; col $1/C #2078b4 target acs')

    run(session,'alias col_t7_doms_1 col $1/A #75c376 target acs; col $1/C #8c6baf target acs; col $1/B #8c96c5 target acs')

    run(session,'alias default_cartoon_and_cols_1  col_t7_doms_1 $1 ; show_this_model $1' )
    
    run(session,'alias default_cartoon_and_cols_2 col_t7_doms_2 $1 ; show_this_model $1' )

    run(session,"alias pic_gamut_lap12 surface; default_cartoon_and_cols_1 #1 ; spin_pic $1_cartoon; surface ; spin_pic $1_surface; default_cartoon_and_cols_1 #1; mlp_bydomain #1 A ; spin_pic $1_A_mlp; default_cartoon_and_cols_1 #1; mlp_bydomain #1 B; spin_pic $1_B_mlp; default_cartoon_and_cols_1 #1; mlp_bydomain #1 C; spin_pic $1_C_mlp; default_cartoon_and_cols_1 #1; col_conservation_bydomain #1 A ; show #1/A surface ; spin_pic $1_A_consurface; default_cartoon_and_cols_1 #1; col_conservation_bydomain #1 B ; show #1/b surface; spin_pic $1_B_consurface; default_cartoon_and_cols_1 #1; col_conservation_bydomain #1 C ; show #1/C surface ; spin_pic $1_C_consurface; default_cartoon_and_cols_1 #1")

    run(session,'alias col_t7_doms_2 col $1/A #75c376 target acs; col $1/B #a7cee2 target acs; col $1/C #2078b4 target acs')

    run(session,"alias pic_gamut_lap34 surface; default_cartoon_and_cols_2 #1 ; spin_pic $1_cartoon; surface ; spin_pic $1_surface; default_cartoon_and_cols_2 #1; mlp_bydomain #1 A ; spin_pic $1_A_mlp; default_cartoon_and_cols_2 #1; mlp_bydomain #1 B; spin_pic $1_B_mlp; default_cartoon_and_cols_2 #1; mlp_bydomain #1 C; spin_pic $1_C_mlp; default_cartoon_and_cols_2 #1; col_conservation_bydomain #1 A ; show #1/A surface ; spin_pic $1_A_consurface; default_cartoon_and_cols_2 #1; col_conservation_bydomain #1 B ; show #1/b surface; spin_pic $1_B_consurface; default_cartoon_and_cols_2 #1; col_conservation_bydomain #1 C ; show #1/C surface ; spin_pic $1_C_consurface; default_cartoon_and_cols_2 #1")

    run(session,"alias multi_pics_full med_licorice ; pic_gamut_lap34 $1_medtrace ; cylinders ; pic_gamut_lap34 $1_cyls")

    run(session,"alias multi_pics_34 thick_licorice ; pic_gamut_lap34 $1_thicktrace")

    run(session,"alias multi_pics_12 thick_licorice ; pic_gamut_lap12 $1_thicktrace")

    return

def register_command(session):
    from chimerax.core.commands import CmdDesc, register
    desc = CmdDesc(synopsis='Start your session off right.')
    register('SetUpSession', desc, SetUpSession, logger=session.logger)
    return


register_command(session)

if __name__ == '__main__':
    from chimerax.core.commands import run
    run(session, 'SetUpSession')