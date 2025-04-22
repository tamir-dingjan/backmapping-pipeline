from enum import Enum


class LipidClass(Enum):
    PC = "PC"
    PE = "PE"
    PI = "PI"
    PA = "PA"
    PS = "PS"
    PG = "PG"
    Cer = "Cer"
    DeoxyCer = "DeoxyCer"
    SM = "SM"
    HexCer = "HexCer"
    Sterol = "Sterol"


unsupported_lipid_classes = {LipidClass.Sterol}
