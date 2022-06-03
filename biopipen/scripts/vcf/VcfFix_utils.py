import re

from biopipen.utils.vcf import *  # noqa: F401, F403


def line_to_obj(line: str):

    for line_obj in (
        HeaderInfo,
        HeaderFormat,
        HeaderFilter,
        HeaderContig,
        HeaderGeneral,
        Fields,
        Variant,
    ):
        if line_obj.is_type(line):
            return line_obj.from_str(line)

    raise ValueError("Unknown line type: {}".format(line))


def handle_obj(obj, fixes: dict):

    for fix in fixes:
        kind = fix.get("kind")
        if kind and obj.kind != kind:
            continue

        id = fix.get("id")
        if id:
            if isinstance(id, str):
                id = [id]

            if obj.kind == "variant" and obj.id not in id:
                continue
            if isinstance(obj, HeaderItem) and obj.get("ID") not in id:
                continue

            return fix["fix"](obj.raw if kind is None else obj)

        regex = fix.get("regex")
        if regex:
            if not re.search(regex, obj.raw):
                continue

            return fix["fix"](obj.raw if kind is None else obj)

        return fix["fix"](obj.raw if kind is None else obj)

    return None


def fix_vcffile(vcffile, outfile, fixes):
    header_append_fixes = []
    variant_append_fixes = []
    modify_fixes = []
    for fix in fixes:
        if fix.get("append") and fix.get("kind") == "variant":
            variant_append_fixes.append(fix)
        elif fix.get("append"):
            header_append_fixes.append(fix)
        else:
            modify_fixes.append(fix)

    with open(vcffile, "r") as fin, open(outfile, "w") as fout:
        for line in fin:
            obj = line_to_obj(line)
            out = handle_obj(obj, modify_fixes)
            if obj.kind == "fields":
                for fix in header_append_fixes:
                    fout.write(str(fix["fix"](None)).rstrip("\n") + "\n")

            if out is False:
                continue
            elif out is None:
                fout.write(str(obj) + "\n")
            else:
                fout.write(str(out).rstrip("\n") + "\n")

        for fix in variant_append_fixes:
            out = fix["fix"](None)
            fout.write(str(out).rstrip("\n") + "\n")
