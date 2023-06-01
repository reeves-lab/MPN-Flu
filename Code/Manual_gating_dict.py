manual_gates_dict = {'B_cells_CD14neg':
                         {'cell': 'B', 'pretty_name': 'B cells'},
                     'Memory_B':
                         {'cell': 'B', 'pretty_name': 'Memory B'},
                     'Plasmablast_38hi27hi':
                         {'cell': 'B', 'pretty_name': 'Plasmablasts'},
                     'Naive_B':
                         {'cell': 'B', 'pretty_name': 'Naive B'},
                     'other_B':
                         {'cell': 'B', 'pretty_name': 'Other B'},
                     'Innate_cells':
                         {'cell': 'I', 'pretty_name': 'Innate cells'},
                     'NK_cells_56hi16hi':
                         {'cell': 'I', 'pretty_name': 'NK 56+16+'},
                     'NK_cell_56hi16neg':
                         {'cell': 'I', 'pretty_name': 'NK 56++16-'},
                     'NK_cells_56hi16lo':
                         {'cell': 'I', 'pretty_name': 'NK 56+16-'},
                     'Intermediate_monocytes':
                         {'cell': 'I', 'pretty_name': 'INT mono'},
                     'Nonclassical_monocytes':
                         {'cell': 'I', 'pretty_name': 'NC mono'},
                     'Classical_monocytes':
                         {'cell': 'I', 'pretty_name': 'CL mono'},
                     'Dendritic_cells':
                         {'cell': 'I', 'pretty_name': 'DC'},
                     'Lin_neg':
                         {'cell': 'I', 'pretty_name': 'Lin-'},
                     'T_cells':
                         {'cell': 'T', 'pretty_name': 'T cells'},
                     'Effector_helper_T':
                         {'cell': 'T', 'pretty_name': 'Effector helper T'},
                     'Naive_helper_T':
                         {'cell': 'T', 'pretty_name': 'Naive helper T'},
                     'CM_helper_T':
                         {'cell': 'T', 'pretty_name': 'CM helper T'},
                     'EM_helper_T':
                         {'cell': 'T', 'pretty_name': 'EM helper T'},
                     'Effector_cytotoxic_T':
                         {'cell': 'T', 'pretty_name': 'Effector cytotoxic T'},
                     'Naive_cytotoxic_T':
                         {'cell': 'T', 'pretty_name': 'Naive cytotoxic T'},
                     'CM_cytotoxic_T':
                         {'cell': 'T', 'pretty_name': 'CM cytotoxic T'},
                     'EM_cytotoxic_T':
                         {'cell': 'T', 'pretty_name': 'EM cytotoxic T'},
                     'DN_T_cells':
                         {'cell': 'T', 'pretty_name': 'DN T'},
                     'NK_T_cells':
                         {'cell': 'T', 'pretty_name': 'NK T'}}

gate_order = {
    'B': {
        'order': ['Naive B', 'Plasmablasts', 'Memory B', 'Other B'],
        'remove': 'B cells'},
    'I': {
        'order': ['CL mono', 'INT mono', 'NC mono', 'DC',
                  'NK (56+16+)', 'NK (56+16-)', 'NK (56++16-)',
                  'Lin-'],
        'remove': 'Innate cells'},
    'T': {
        'order': ['Naive helper T', 'Effector helper T', 'CM helper T', 'EM helper T',
                  'Naive cytotoxic T', 'Effector cytotoxic T', 'CM cytotoxic T', 'EM cytotoxic T',
                  'DN T', 'NK T'],
        'remove': 'T cells'}
}
