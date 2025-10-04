# meal_satiety.py
import re, sys, time
from typing import Dict, Any, Optional, List
import requests

# ---------- Utils nutriments ----------
def _get_kcal(n: Dict[str, Any]) -> Optional[float]:
    if not n: return None
    kcal = n.get("energy-kcal_100g")
    if kcal is not None:
        try: return float(kcal)
        except: return None
    kj = n.get("energy_100g")
    if kj is not None:
        try: return float(kj) / 4.184
        except: return None
    return None

def _get(n: Dict[str, Any], key: str, alt: Optional[str] = None, default: float = 0.0) -> float:
    if not n: return default
    v = n.get(key)
    if v is None and alt: v = n.get(alt)
    try: return float(v) if v is not None else default
    except: return default

def _salt_from_nutriments(n: Dict[str, Any]) -> float:
    if not n: return 0.0
    if "salt_100g" in n:
        try: return float(n["salt_100g"])
        except: pass
    if "sodium_100g" in n:
        try: return float(n["sodium_100g"]) * 2.5
        except: pass
    return 0.0

def fetch_off(barcode: str) -> Dict[str, Any]:
    url = f"https://world.openfoodfacts.org/api/v0/product/{barcode}.json"
    r = requests.get(url, timeout=15)
    r.raise_for_status()
    js = r.json()
    if js.get("status") != 1:
        raise ValueError("Produit introuvable sur OpenFoodFacts")
    return js["product"]

# ---------- Heuristiques ----------
def detect_is_liquid(product: Dict[str, Any]) -> str:
    tags = set(product.get("categories_tags", []) or [])
    if any(t in tags for t in ["en:beverages","en:carbonated-drinks","en:juices",
                               "en:plant-based-milks","en:teas","en:coffees","en:waters"]):
        return "liquid"
    if any(t in tags for t in ["en:dairy-desserts","en:puddings","en:compotes","en:spreads"]):
        return "semisolid"
    return "solid"

def detect_upf_penalty(product: Dict[str, Any]) -> bool:
    nova = product.get("nova_group")
    additives = product.get("additives_tags") or []
    try: nova = int(nova) if nova is not None else None
    except: nova = None
    return (nova == 4) or (len(additives) >= 4)

# ---------- IS (g²/kcal) ----------
def compute_is_from_nutriments(n: Dict[str, Any], shape: str, upf_pen: bool) -> float:
    E   = _get_kcal(n)
    if not E or E <= 0: raise ValueError("Calories manquantes")
    P   = _get(n, "proteins_100g")
    F   = _get(n, "fiber_100g")
    L   = _get(n, "fat_100g")
    SFA = _get(n, "saturated-fat_100g", "saturated_fat_100g")
    SUt = _get(n, "sugars_100g")
    SA  = _get(n, "added-sugars_100g")
    if SA == 0.0: SA = _get(n, "sugars-added_100g")
    SUG = SA if SA > 0 else SUt
    SALT = _salt_from_nutriments(n)

    P_c=min(P,20.0); F_c=min(F,12.0); UF=max(L-SFA,0.0)
    SUG_t = 15 + 0.5*(SUG-15) if SUG>15 else SUG
    numerator = (1.0*P_c)+(0.8*F_c)+(0.5*UF)-(1.1*SUG_t)-(0.4*SFA)-(0.2*SALT)
    IS = numerator / (E/100.0)
    mod = 1.0
    if shape == "liquid": mod *= 0.85
    elif shape == "semisolid": mod *= 0.95
    if upf_pen: mod *= 0.90
    return IS*mod

def compute_satiety_index_v2(product: Dict[str, Any]) -> Dict[str, Any]:
    n = product.get("nutriments", {}) or {}
    E = _get_kcal(n); P=_get(n,"proteins_100g"); F=_get(n,"fiber_100g"); L=_get(n,"fat_100g")
    SFA=_get(n,"saturated-fat_100g","saturated_fat_100g"); SUt=_get(n,"sugars_100g")
    SA=_get(n,"added-sugars_100g"); SA = SA if SA>0 else _get(n,"sugars-added_100g")
    SALT = _salt_from_nutriments(n); CARB=_get(n,"carbohydrates_100g")

    shape = detect_is_liquid(product); upf = detect_upf_penalty(product)
    ISm = compute_is_from_nutriments(n, shape, upf)

    fields = sum([1 if E else 0, 1 if P else 0, 1 if F else 0, 1 if L else 0,
                  1 if SFA else 0, 1 if (SA or SUt) else 0, 1 if SALT else 0])
    conf = round(fields/7.0,2)

    return {"E_kcal_100g":round(E or 0,1),"P":P,"F":F,"L":L,"SFA":SFA,
            "SUG_total":SUt,"SUG_added":(SA if SA else None),"CARB":CARB,
            "SALT":SALT,"shape":shape,"nova_group":product.get("nova_group"),
            "additives_n":len(product.get("additives_tags") or []),
            "I_v2_mod_g2kcal":ISm,"confidence":conf,
            "product_name":product.get("product_name"),"brands":product.get("brands")}

def is_band(v: float) -> str:
    return "Excellent" if v>6 else "Bon" if v>3 else "Moyen" if v>1 else "Faible" if v>0 else "Anti-satiété"
def se_band(v: float) -> str:
    return "Repas très rassasiant" if v>20 else "Bon effet coupe-faim" if v>8 else "Satiété modérée" if v>3 else "Faible satiété" if v>1 else "Négligeable (plaisir)"

# ---------- Estimation faim ----------
def estimate_time_to_hunger(se_cal: float, kcal_portion: float, shape: str, confidence: float):
    if se_cal < 1: base = 0.8
    elif se_cal < 3: base = 1.5 + 0.3*(se_cal-1)
    elif se_cal < 8: base = 2.1 + 0.3*(se_cal-3)
    elif se_cal < 20: base = 3.6 + 0.12*(se_cal-8)
    else: base = 5.0 + 0.08*(min(se_cal,35)-20)
    base *= min(max(kcal_portion/600.0,0.6),1.2)
    if shape == "liquid": base *= 0.85
    elif shape == "semisolid": base *= 0.95
    base = max(0.5, min(base, 7.0))
    delta = max(0.5, (1.2 - confidence) * 1.25)
    return base, max(0.25, base-delta), min(8.0, base+delta)

# ---------- Portion ----------
DEFAULT_PORTIONS = {"repas":300.0,"snack":50.0,"plaisir":15.0}
def ask_portion_interactive(default_g: float=100.0) -> float:
    prompt=f"Portion : "
    while True:
        ans=input(prompt).strip().lower()
        if ans=="": return float(default_g)
        if ans in DEFAULT_PORTIONS: return DEFAULT_PORTIONS[ans]
        m=re.match(r"^\s*(\d+(\.\d+)?)\s*g?\s*$",ans)
        if m:
            v=float(m.group(1))
            if v>0: return v
        print("Entrée invalide. Ex: 250 ou repas|snack|plaisir")

# ---------- Saisie manuelle ----------
def manual_entry_to_product() -> (Dict[str,Any], float):
    print("\n=== Saisie MANUELLE par PORTION ===")
    portion=float(input("Poids portion (g)       : ").strip())
    kcal   =float(input("Énergie (kcal)          : ").strip())
    P      =float(input("Protéines (g)           : ").strip())
    F      =float(input("Fibres (g)              : ").strip())
    L      =float(input("Lipides (g)             : ").strip())
    SFA    =float(input("Saturés (g)             : ").strip())
    SUG    =float(input("Sucres totaux (g)       : ").strip())
    SALT   =float(input("Sel (g)                 : ").strip())
    shape  =(input("Forme (solid/semisolid/liquid) [solid] : ").strip().lower() or "solid")
    fac100 = 100.0/max(1.0,portion)
    n={"energy-kcal_100g":kcal*fac100,"proteins_100g":P*fac100,"fiber_100g":F*fac100,
       "fat_100g":L*fac100,"saturated-fat_100g":SFA*fac100,"sugars_100g":SUG*fac100,"salt_100g":SALT*fac100}
    prod={"product_name":"Entrée manuelle","brands":"","nutriments":n,"categories_tags":[], "nova_group":None,"additives_tags":[]}
    if shape=="liquid": prod["categories_tags"]=["en:beverages"]
    elif shape=="semisolid": prod["categories_tags"]=["en:puddings"]
    return prod, portion

# ---------- Webcam (1 scan puis retour) ----------
def webcam_scan_once() -> Optional[str]:
    try:
        import cv2
        from pyzbar.pyzbar import decode
    except Exception:
        print("Webcam indisponible (opencv/pyzbar manquants)")
        return None
    cap=cv2.VideoCapture(0)
    if not cap.isOpened():
        print("Impossible d'ouvrir la webcam"); return None
    print("Visez un code-barres… (ESC/q pour annuler)")
    last=None; last_t=0
    try:
        while True:
            ok, frame = cap.read()
            if not ok: break
            for bc in decode(frame):
                code = bc.data.decode("utf-8")
                if code!=last or (time.time()-last_t)>1.2:
                    last, last_t = code, time.time()
                    cv2.imshow("Scan", frame)
                    cap.release(); cv2.destroyAllWindows()
                    return code
            cv2.imshow("Scan", frame)
            k = cv2.waitKey(1) & 0xFF
            if k in (27, ord('q')): break
    finally:
        try:
            cap.release(); cv2.destroyAllWindows()
        except: pass
    return None

# ---------- Qualité (commentaire bref) ----------
def quality_comment(is_value: float, res: dict) -> str:
    P=res.get("P") or 0.0; F=res.get("F") or 0.0; SFA=res.get("SFA") or 0.0
    SALT=res.get("SALT") or 0.0; SUG=res.get("SUG_total") or 0.0; E=res.get("E_kcal_100g") or 0.0
    lowP=P<8; lowF=F<3; hiSFA=SFA>5; hiSALT=SALT>1.2; hiSUG=SUG>10; dense=E>300
    if is_value>6: head="Qualité excellente."
    elif is_value>3: head="Qualité bonne."
    elif is_value>1: head="Qualité moyenne."
    elif is_value>0: head="Qualité faible."
    else: head="Qualité très faible."
    tips=[]
    if lowP: tips.append("↑ protéines");
    if lowF: tips.append("↑ fibres")
    if hiSFA: tips.append("↓ saturés")
    if hiSALT: tips.append("↓ sel")
    if hiSUG: tips.append("↓ sucres")
    if dense: tips.append("↓ densité énergétique")
    return head + (" " + " ; ".join(tips) if tips else "")

# ---------- Run pour un item ----------
def run_on_product(product: Dict[str,Any], portion_g: float, targets=(25.0,8.0,600.0)):
    res = compute_satiety_index_v2(product)
    Fq = portion_g/100.0
    IS = res["I_v2_mod_g2kcal"]
    SE = IS * Fq
    TP,TF,TK = targets
    Pp=(res["P"] or 0)*Fq; Fp=(res["F"] or 0)*Fq; Kp=res["E_kcal_100g"]*Fq
    prot_factor=min(Pp/TP,1.0) if TP>0 else 1.0
    fib_factor =min(Fp/TF,1.0) if TF>0 else 1.0
    kcal_factor=min(Kp/TK,1.0) if TK>0 else 1.0
    bread_penalty=1.0
    name=(product.get("product_name") or "").lower()
    cats=" ".join(product.get("categories_tags") or [])
    if ("sandwich" in name or "en:meals" in cats or "en:prepared-meals" in cats) and ((res["F"] or 0)<3.0):
        bread_penalty=0.85
    SE_cal = SE * prot_factor * fib_factor * kcal_factor * bread_penalty

    eta, lo, hi = estimate_time_to_hunger(SE_cal, Kp, res["shape"], res["confidence"])

    # apports pour la portion
    item = {
        "name": f"{product.get('product_name')} [{product.get('brands') or ''}]".strip(),
        "portion_g": portion_g,
        "kcal": Kp,
        "P": Pp, "F": Fp,
        "CARB": (res.get("CARB") or 0)*Fq,
        "SUG": (res.get("SUG_total") or 0)*Fq,
        "L": (res.get("L") or 0)*Fq,
        "SFA": (res.get("SFA") or 0)*Fq,
        "SALT": (res.get("SALT") or 0)*Fq,
        "IS": IS, "SE_cal": SE_cal, "shape": res["shape"],
        "conf": res["confidence"], "quality": quality_comment(IS, res)
    }

    # affichage court par item
    print("\n— Item —")
    print(f"{item['name']}")
    print(f"Portion {portion_g:.0f} g | {item['kcal']:.0f} kcal | IS {IS:.2f} | SE_cal {SE_cal:.2f} ({se_band(SE_cal)})")
    print(f"Prot {item['P']:.1f} g • Fib {item['F']:.1f} g • Glu {item['CARB']:.1f} g (Sucres {item['SUG']:.1f} g)")
    print(f"Lip {item['L']:.1f} g (Sat {item['SFA']:.1f} g) • Sel {item['SALT']:.2f} g")
    print(f"Qualité: {item['quality']}")

    return item

# ---------- Agrégation repas ----------
def aggregate_meal(items: List[Dict[str,Any]]):
    if not items: return
    tot_g = sum(i["portion_g"] for i in items)
    tot_kcal = sum(i["kcal"] for i in items)
    tot_P = sum(i["P"] for i in items)
    tot_F = sum(i["F"] for i in items)
    tot_C = sum(i["CARB"] for i in items)
    tot_SUG = sum(i["SUG"] for i in items)
    tot_L = sum(i["L"] for i in items)
    tot_SFA = sum(i["SFA"] for i in items)
    tot_SALT = sum(i["SALT"] for i in items)

    # profil /100 g du repas (pondération au poids)
    if tot_g <= 0: return
    fac100 = 100.0/tot_g
    n = {
        "energy-kcal_100g": tot_kcal * fac100,
        "proteins_100g":    tot_P    * fac100,
        "fiber_100g":       tot_F    * fac100,
        "fat_100g":         tot_L    * fac100,
        "saturated-fat_100g": tot_SFA* fac100,
        "sugars_100g":      tot_SUG  * fac100,
        "salt_100g":        tot_SALT * fac100,
        "carbohydrates_100g": tot_C  * fac100
    }
    # forme : on considère solide, sauf si majorité de liquides
    liquid_ratio = sum(1 for i in items if i["shape"]=="liquid")/len(items)
    shape = "liquid" if liquid_ratio>0.6 else "solid"
    IS_meal = compute_is_from_nutriments(n, shape, upf_pen=False)
    Fq = tot_g/100.0
    # calibration repas (mêmes cibles)
    TP,TF,TK = 25.0, 8.0, 600.0
    prot_factor = min(tot_P/TP,1.0)
    fib_factor  = min(tot_F/TF,1.0)
    kcal_factor = min(tot_kcal/TK,1.0)
    SE_cal_meal = IS_meal * Fq * prot_factor * fib_factor * kcal_factor
    eta, lo, hi = estimate_time_to_hunger(SE_cal_meal, tot_kcal, shape, confidence=1.0)

    print("\n==============================")
    print("RÉCAP REPAS (tous les items)")
    print(f"Total: {tot_g:.0f} g • {tot_kcal:.0f} kcal")
    print(f"Prot {tot_P:.1f} g • Fib {tot_F:.1f} g • Glu {tot_C:.1f} g (Sucres {tot_SUG:.1f} g)")
    print(f"Lip {tot_L:.1f} g (Sat {tot_SFA:.1f} g) • Sel {tot_SALT:.2f} g")
    print(f"IS du repas: {IS_meal:.2f} g²/kcal ({is_band(IS_meal)})")
    print(f"SE_cal du repas: {SE_cal_meal:.2f} → {se_band(SE_cal_meal)}")
    print(f"≈ faim dans ~{eta:.1f} h (intervalle {lo:.1f}–{hi:.1f} h)")

# ---------- Main ----------
def main():
    try:
        n = int(input("Combien d'aliments veux-tu évaluer pour ce repas ? ").strip())
    except:
        n = 1

    items = []
    for idx in range(1, n+1):
        print(f"\n=== Item {idx}/{n} ===")
        print("1) Code-barres  2) Saisie manuelle  3) Webcam")
        choice = (input("Choix [1/2/3] : ").strip() or "1")
        try:
            if choice == "2":
                product, portion = manual_entry_to_product()
            elif choice == "3":
                code = webcam_scan_once()
                if not code:
                    print("Scan annulé. Passage en saisie manuelle.")
                    product, portion = manual_entry_to_product()
                else:
                    product = fetch_off(code)
                    portion = ask_portion_interactive(100.0)
            else:
                code = input("Code-barres : ").strip()
                product = fetch_off(code)
                portion = ask_portion_interactive(100.0)

            item = run_on_product(product, portion)
            items.append(item)

        except Exception as e:
            print(f"[Item {idx}] Erreur: {e}")

    aggregate_meal(items)

if __name__ == "__main__":
    main()
