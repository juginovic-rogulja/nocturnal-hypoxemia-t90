"""
Diagnostic code definitions for every outcome, v8 (2026-09-12). Lane X2, steps 040/041.

Derived from numbers/disease_definitions.py (v7, left untouched on disk) with the six list fixes of decision 5 and the
four added outcomes of decision 11 (DECISIONS_2026-09-12.md; RERUN_PLAN_V8.md item D3). Every other entry is copied
verbatim; test_definitions_v8.py proves that against the v7 module at run time.

Matching rule (unchanged from v7): strip the decimal point from condition_source_value, upper-case, prefix match.
"I50" matches I50, I50.9, I5031. New in v8: an optional per-condition EXCLUDE list of SQL LIKE patterns applied after
the include match ("_" is any single character), used by type 2 diabetes to drop the ICD-9 type 1 fifth digits.

Code meanings are written next to every changed or new prefix.

Decision 5 fixes
  fatty liver     571.8 other chronic nonalcoholic liver disease, 571.9 unspecified chronic liver disease without alcohol,
                  K76.0 fatty (change of) liver not elsewhere classified, K75.8 other specified inflammatory liver diseases
                  (K75.81 nonalcoholic steatohepatitis), K75.9 inflammatory liver disease unspecified.
                  REMOVED 571.5 cirrhosis of liver without mention of alcohol and K74.6 other and unspecified cirrhosis of liver:
                  both are cirrhosis and were also in the cirrhosis list, so a cirrhosis patient counted as fatty liver too.
  nocturia        788.43 nocturia (was 788.19 other dysuria, which has zero rows in the cache), R35.1 nocturia.
  dementia        adds 331.0 Alzheimer's disease (ICD-9 twin of G30, which was already in).
  type 2 diabetes 250 diabetes mellitus MINUS the type 1 fifth digits 250.x1 (type I, not stated as uncontrolled) and
                  250.x3 (type I, uncontrolled); E11 type 2 diabetes mellitus, E13 other specified diabetes mellitus.
                  E10 type 1 diabetes was never in the list.
  peripheral artery disease  440.2 atherosclerosis of native arteries of the extremities, 443.9 peripheral vascular disease
                  unspecified, 443.81 peripheral angiopathy in diseases classified elsewhere, I70 atherosclerosis,
                  I73.9 peripheral vascular disease unspecified. REMOVED 443.2 other arterial dissection (not atherosclerotic).
  M80             osteoporosis with current pathological fracture: stays in osteoporosis, leaves fracture (M80 was in both, so one
                  diagnosis produced two outcome events). ALTERNATIVE (not applied): M80 into fracture and out of osteoporosis;
                  check_outcomes_v8.py prints the incident counts under both placements.

Decision 11 additions
  polycythemia    238.4 polycythemia vera, 289.0 secondary polycythemia; D45 polycythemia vera, D75.1 secondary polycythemia.
  parkinson       332.0 paralysis agitans (Parkinson's disease); G20 Parkinson's disease.
                  EXCLUDED by omission: 332.1 secondary parkinsonism and G21 secondary parkinsonism (drug-induced, vascular, other).
  ild             515 postinflammatory pulmonary fibrosis; 516 other alveolar and parietoalveolar pneumonopathy (516.3x idiopathic
                  interstitial pneumonias with 516.31 idiopathic pulmonary fibrosis, 516.8 other specified, 516.9 unspecified);
                  J84 other interstitial pulmonary diseases (J84.1x incl. J84.112 idiopathic pulmonary fibrosis, J84.9 unspecified).
  vent_arrhythmia_arrest  427.1 paroxysmal ventricular tachycardia, 427.41 ventricular fibrillation, 427.42 ventricular flutter
                  (prefix 4274 covers both), 427.5 cardiac arrest; I47.2 ventricular tachycardia, I49.0 ventricular fibrillation and
                  flutter, I46 cardiac arrest (I46.2 due to underlying cardiac condition, I46.8 other underlying, I46.9 unspecified).

Decision of 14 Sept 2026 (v8.1), negative-control swap
  cataract -> hernia_ing (550; K40 inguinal hernia) and back_pain -> alopecia (704.0; L63, L64, L65). The two that left sat above 1
  per SD on the v8 tables (1.067 and 1.074) and each has a pathway to the exposure (diabetes, obesity). The two that entered were
  screened in August as mechanical or hereditary conditions with no described pathway from nocturnal oxygenation. Straight swap:
  same roles (control and ranked), so 57 conditions, 52 ranked, 54 diseases and death are unchanged. SWAPPED_CONTROLS_2026_09_14.
"""

# key: (label, negative_control, [icd9 prefixes], [icd10 prefixes])
DISEASES = {
    # ---------------------------------------------------------------- cardiovascular
    "hf":            ("Heart failure", False, ["428"], ["I50", "I110", "I130", "I132"]),
    "ihd":           ("Ischemic heart disease", False, ["410", "411", "412", "413", "414"],
                      ["I20", "I21", "I22", "I23", "I24", "I25"]),
    "mi":            ("Myocardial infarction", False, ["410", "412"], ["I21", "I22", "I252"]),
    "afib":          ("Atrial fibrillation", False, ["42731", "42732"], ["I48"]),
    "stroke_any":    ("Stroke", False, ["430", "431", "432", "433", "434", "436"],
                      ["I60", "I61", "I62", "I63", "I64"]),
    # v8: 443.2 (other arterial dissection) removed
    "pad":           ("Peripheral artery disease", False, ["4402", "4439", "44381"],
                      ["I70", "I739"]),
    "htn2":          ("Hypertension", False, ["401", "402", "403", "404", "405"],
                      ["I10", "I11", "I12", "I13", "I15"]),
    "pulm_htn":      ("Pulmonary hypertension", False, ["4160", "4168", "4169"],
                      ["I270", "I272", "I2720", "I2721", "I2789", "I279"]),
    "vte":           ("Venous thromboembolism", False, ["4151", "4534", "4538", "4539"],
                      ["I26", "I82"]),
    # v8 NEW: ventricular tachycardia, ventricular fibrillation or flutter, cardiac arrest
    "vent_arrhythmia_arrest": ("Ventricular arrhythmia or cardiac arrest", False,
                      ["4271", "4274", "4275"], ["I472", "I490", "I46"]),
    # ---------------------------------------------------------------- respiratory
    "copd2":         ("COPD", False, ["491", "492", "496"], ["J41", "J42", "J43", "J44"]),
    "asthma":        ("Asthma", False, ["493"], ["J45"]),
    "resp_failure":  ("Respiratory failure", False, ["51881", "51882", "51883", "51884", "7991"],
                      ["J96", "R092"]),
    "pneumonia":     ("Pneumonia", False, ["480", "481", "482", "483", "485", "486", "507"],
                      ["J12", "J13", "J14", "J15", "J16", "J17", "J18", "J69"]),
    "obesity_hypovent": ("Obesity hypoventilation", False, ["27803"], ["E662"]),
    # v8 NEW: interstitial lung disease
    "ild":           ("Interstitial lung disease", False, ["515", "516"], ["J84"]),
    # ---------------------------------------------------------------- metabolic
    # v8: ICD-9 type 1 fifth digits excluded through EXCLUDE["diabetes"] below
    "diabetes":      ("Type 2 diabetes", False, ["250"], ["E11", "E13"]),
    "obesity":       ("Obesity", False, ["2780"], ["E66"]),
    "dyslipid":      ("Dyslipidemia", False, ["272"], ["E78"]),
    "gout":          ("Gout", False, ["274"], ["M10", "M1A"]),
    "thyroid_dis":   ("Thyroid disease", False, ["240", "241", "242", "243", "244", "245", "246"],
                      ["E00", "E01", "E02", "E03", "E04", "E05", "E06", "E07"]),
    # ---------------------------------------------------------------- renal, hepatic
    "ckd":           ("Chronic kidney disease", False, ["585"], ["N18"]),
    "aki":           ("Acute kidney injury", False, ["584"], ["N17"]),
    # v8: 571.5 and K74.6 (cirrhosis) removed
    "nafld":         ("Fatty liver disease", False, ["5718", "5719"], ["K760", "K758", "K759"]),
    "cirrhosis":     ("Cirrhosis", False, ["5712", "5715", "5716"], ["K703", "K717", "K74"]),
    # ---------------------------------------------------------------- infection
    "sepsis":        ("Sepsis", False, ["038", "78552", "99591", "99592"],
                      ["A40", "A41", "R652", "R6520", "R6521"]),
    "cellulitis":    ("Cellulitis", False, ["681", "682"], ["L03"]),
    "uti":           ("Urinary tract infection", False, ["5990"], ["N390"]),
    # ---------------------------------------------------------------- neurologic, psychiatric
    # v8: 331.0 Alzheimer's disease added
    "dementia":      ("Dementia", False, ["290", "2941", "3310", "3312"],
                      ["F00", "F01", "F02", "F03", "G30", "G31"]),
    # v8 NEW: Parkinson's disease (secondary parkinsonism 332.1 / G21 excluded)
    "parkinson":     ("Parkinson's disease", False, ["3320"], ["G20"]),
    "epilepsy":      ("Epilepsy", False, ["345"], ["G40"]),
    "migraine2":     ("Migraine", False, ["346"], ["G43"]),
    "depression":    ("Depression", False, ["2962", "2963", "3004", "311"],
                      ["F32", "F33", "F341"]),
    "anxiety2":      ("Anxiety", False, ["3000", "3002", "30002", "30023"], ["F40", "F41"]),
    "bipolar":       ("Bipolar disorder", False, ["2960", "2964", "2965", "2966", "2967"],
                      ["F30", "F31"]),
    "adhd":          ("ADHD", False, ["314"], ["F90"]),
    "insomnia":      ("Insomnia", False, ["30742", "78052"], ["F510", "G470"]),
    "rls_plmd":      ("Restless legs syndrome", False, ["33394", "78058"], ["G2581", "G4761"]),
    # ---------------------------------------------------------------- musculoskeletal
    "osteoarthritis": ("Osteoarthritis", False, ["715"], ["M15", "M16", "M17", "M18", "M19"]),
    # v8.1 (14 Sept 2026): back pain left the control panel (its per-SD ratio sat above 1 with the interval excluding 1, and back
    # pain tracks obesity and mechanical load). Alopecia takes its slot, as a negative control AND as a ranked outcome, so every
    # count of the paper is unchanged. Code meanings: 704.0 alopecia; L63 alopecia areata, L64 androgenic alopecia, L65 other
    # nonscarring hair loss (the August screen list, predominantly androgenetic, hormonal and hereditary).
    "alopecia":      ("Alopecia", True, ["7040"], ["L63", "L64", "L65"]),
    # M80 stays here (default placement, decision 5)
    "osteoporosis":  ("Osteoporosis", False, ["7330"], ["M80", "M81"]),
    # v8: M80 removed (it stays in osteoporosis)
    "fracture":      ("Fracture", False, ["800", "801", "802", "805", "807", "808", "810",
                                          "812", "813", "820", "821", "823", "824"],
                      ["S02", "S12", "S22", "S32", "S42", "S52", "S62", "S72", "S82", "S92"]),
    "fibromyalgia":  ("Fibromyalgia", False, ["7291"], ["M797"]),
    "falls":         ("Falls", False, ["E8880", "E8881", "E8888", "E8889", "E8840", "E8841",
                                   "E8842", "E8843", "E8844", "E8845", "E8846", "E8849"],
                  ["W01", "W06", "W07", "W08", "W10",
                                                         "W18", "W19", "R296"]),
    # ---------------------------------------------------------------- sensory
    # v8.1 (14 Sept 2026): cataract left the control panel (its per-SD ratio crossed 1, and cataract tracks diabetes and age).
    # Inguinal hernia takes its slot, as a negative control AND as a ranked outcome. Code meanings: 550 inguinal hernia (all
    # fifth digits); K40 inguinal hernia (the August screen list, a mechanical defect of the inguinal canal).
    "hernia_ing":    ("Inguinal hernia", True, ["550"], ["K40"]),
    "glaucoma":      ("Glaucoma", True, ["365"], ["H40", "H42"]),
    "hearing_loss2": ("Hearing loss", False, ["389"], ["H90", "H91"]),
    # ------------------------------------------------- replacement negative controls, 2026-08-07
    "dermatitis_contact": ("Contact dermatitis", True, ["692"], ["L23", "L24", "L25"]),
    "haemorrhoids":  ("Hemorrhoids", True, ["455"], ["I84", "K64"]),
    # ---------------------------------------------------------------- other
    "anaemia":       ("Anemia", False, ["280", "281", "285"], ["D50", "D51", "D52", "D53",
                                                               "D63", "D64"]),
    # v8 NEW: polycythemia vera and secondary polycythemia
    "polycythemia":  ("Polycythemia", False, ["2384", "2890"], ["D45", "D751"]),
    "gerd":          ("Gastroesophageal reflux", False, ["53081"], ["K21"]),
    "cancer_any":    ("Cancer", False, ["14", "15", "16", "17", "18", "19", "200", "201",
                                        "202", "203", "204", "205", "206", "207", "208"],
                      ["C"]),
    "prostate_ca":   ("Prostate cancer", False, ["185"], ["C61"]),
    "erectile":      ("Erectile dysfunction", False, ["60784"], ["N52"]),
    # v8: 788.43 nocturia (was 788.19 other dysuria)
    "nocturia":      ("Nocturia", False, ["78843"], ["R351"]),
    "cvd":           ("Cardiovascular composite", False, [], []),   # composite, built below
}

# Optional per-condition exclusion patterns (SQL LIKE on the dot-stripped upper-cased code; "_" = one character),
# applied AFTER the include prefixes. Type 2 diabetes: drop ICD-9 250.x1 and 250.x3, the type 1 fifth digits.
EXCLUDE = {
    "diabetes": ["250_1", "250_3"],
}

# the composite is the union of its components rather than a code list of its own (unchanged; PAD's fix flows in)
CVD_COMPONENTS = ["hf", "ihd", "mi", "stroke_any", "pad"]

# what v8 changed (decision 5) and added (decision 11); the disclosure sentence reads ADDED_2026_09_12
CHANGED_2026_09_12 = {"nafld", "nocturia", "dementia", "diabetes", "pad", "fracture"}
ADDED_2026_09_12 = {"polycythemia", "parkinson", "ild", "vent_arrhythmia_arrest"}
# conditions the v8 builder recomputes from the cache: the changed six, the four new, osteoporosis (same list, the M80
# twin of fracture, recomputed as an identity control) and the composite (PAD is a component)
# v8.1 (14 Sept 2026, Alen): the two controls whose per-SD ratio sat above 1 were swapped for two from the August screen,
# each taking the old one's roles (negative control and ranked outcome); disclosed in Methods with the decision date.
SWAPPED_CONTROLS_2026_09_14 = {"cataract": "hernia_ing", "back_pain": "alopecia"}
SWAPPED_IN_2026_09_14 = set(SWAPPED_CONTROLS_2026_09_14.values())
RECOMPUTE_SET = CHANGED_2026_09_12 | ADDED_2026_09_12 | {"osteoporosis", "cvd"} | SWAPPED_IN_2026_09_14
M80_ALTERNATIVE = "M80 into fracture and out of osteoporosis (not applied; counts under both placements in the report)"

_DROPPED_CONTROLS = {"thyroid_dis", "fracture", "osteoarthritis"}
assert all(DISEASES[k][1] is False for k in _DROPPED_CONTROLS), _DROPPED_CONTROLS
NEGATIVE_CONTROLS = [k for k, v in DISEASES.items() if v[1] and k not in _DROPPED_CONTROLS]
assert NEGATIVE_CONTROLS == ["alopecia", "hernia_ing", "glaucoma",
                             "dermatitis_contact", "haemorrhoids"], NEGATIVE_CONTROLS   # v8.1 panel
NEGATIVE_CONTROL_LABELS = [DISEASES[k][0] for k in NEGATIVE_CONTROLS]

# v8 F2 fix 2026-09-12 (integrity gate 6): no typed floor. The floor is a RESULT, written to numbers/negcontrols_final.json
# by make_negcontrols_final.py (chain step 111b, from bdsp_diseases_v3.csv of step 111) and read from there, so the
# scripts that import this constant (q2_noapnea_freq/model.py step 133, habitual_outcomes_v2/attack.py step 151) see
# the same floor as the json readers (splitstyle.py, model_habitual.py). The earlier comment ("refreshed by
# run_primary_unpenalized.py") was wrong: that script never touched this constant, so the v7 value 1.063 sat here.
# None only on a fresh tree where the writer itself imports this module before the json exists.
try:
    with open(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)),
                                         "negcontrols_final.json")) as _fh:
        CONFOUNDING_FLOOR_PER_SD = float(__import__("json").load(_fh)["confounding_floor_per_sd"])
except FileNotFoundError:
    CONFOUNDING_FLOOR_PER_SD = None

RANKING_EXCLUDE = {"dermatitis_contact", "haemorrhoids"}   # unchanged (plan, step 040)

ORGAN_GROUP = {
    "hf": "Cardiac", "ihd": "Cardiac", "mi": "Cardiac", "afib": "Cardiac",
    "pulm_htn": "Cardiac", "htn2": "Cardiac", "pad": "Cardiac", "vte": "Cardiac",
    "vent_arrhythmia_arrest": "Cardiac",
    "stroke_any": "Neuro/psych",
    "copd2": "Respiratory", "asthma": "Respiratory", "resp_failure": "Respiratory",
    "pneumonia": "Respiratory", "obesity_hypovent": "Respiratory", "ild": "Respiratory",
    "diabetes": "Metabolic", "obesity": "Metabolic", "dyslipid": "Metabolic",
    "gout": "Metabolic", "thyroid_dis": "Metabolic",
    "ckd": "Kidney", "aki": "Kidney",
    "nafld": "Liver", "cirrhosis": "Liver",
    "sepsis": "Infection", "cellulitis": "Infection", "uti": "Infection",
    "dementia": "Neuro/psych", "parkinson": "Neuro/psych", "epilepsy": "Neuro/psych", "migraine2": "Neuro/psych",
    "depression": "Neuro/psych", "anxiety2": "Neuro/psych", "bipolar": "Neuro/psych",
    "adhd": "Neuro/psych", "insomnia": "Neuro/psych", "rls_plmd": "Neuro/psych",
    "osteoarthritis": "Sensory/MSK",
    "osteoporosis": "Sensory/MSK", "fracture": "Sensory/MSK",
    "fibromyalgia": "Sensory/MSK", "falls": "Sensory/MSK",
    "glaucoma": "Sensory/MSK", "hearing_loss2": "Sensory/MSK",
    "anaemia": "Other", "polycythemia": "Other", "gerd": "Other", "cancer_any": "Other", "prostate_ca": "Other",
    "erectile": "Other", "nocturia": "Other",
    "dermatitis_contact": "Other", "haemorrhoids": "Other", "hernia_ing": "Other", "alopecia": "Other",
    "cvd": "Composite",
}
assert set(ORGAN_GROUP) == set(DISEASES), set(DISEASES) ^ set(ORGAN_GROUP)

CIRCULAR = ["insomnia", "rls_plmd", "nocturia", "epilepsy"]   # unchanged
