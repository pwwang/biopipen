from __future__ import annotations
from typing import Any, Sequence


class HeaderItem(dict):
    """The base class of header items"""
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.raw = None

    @classmethod
    def from_str(cls, line: str):
        obj = cls()
        obj.raw = line

        line = line.rstrip("\r\n")
        items = line[line.find("<") + 1 : -1].split(",", 3)
        for item in items:
            key, value = item.split("=", 1)
            if key == "Description":
                value = value[1:-1]
            obj[key] = value

        return obj

    def __setattr__(self, name: str, value: Any) -> None:
        return super().__setitem__(name, value)

    def __getattr__(self, name: str) -> Any:
        return super().__getitem__(name)


class HeaderInfo(HeaderItem):
    """The INFO items in the header"""

    kind = "info"

    def __str__(self):
        return (
            f"##INFO=<ID={self['ID']},"
            f"Number={self['Number']},"
            f"Type={self['Type']},"
            f"Description=\"{self['Description']}\">"
        )

    @staticmethod
    def is_type(raw: str) -> bool:
        return raw.startswith("##INFO")


class HeaderFormat(HeaderItem):
    """The FORMAT items in the header"""

    kind = "format"

    def __str__(self):
        return (
            f"##FORMAT=<ID={self['ID']},"
            f"Number={self['Number']},"
            f"Type={self['Type']},"
            f"Description=\"{self['Description']}\">"
        )

    @staticmethod
    def is_type(raw: str) -> bool:
        return raw.startswith("##FORMAT")


class HeaderFilter(HeaderItem):
    """The FILTER items in the header"""

    kind = "filter"

    def __str__(self):
        return (
            f"##FILTER=<ID={self['ID']},"
            f"Description=\"{self['Description']}\">"
        )

    @staticmethod
    def is_type(raw: str) -> bool:
        return raw.startswith("##FILTER")


class HeaderContig(HeaderItem):
    """The contig items in the header"""

    kind = "contig"

    def __str__(self):
        return f"##contig=<ID={self['ID']}," f"length={self['length']}>"

    @staticmethod
    def is_type(raw: str) -> bool:
        return raw.startswith("##contig")


class HeaderGeneral(HeaderItem):
    """The general items in the header"""

    kind = "header"

    @classmethod
    def from_str(cls, line: str):
        obj = cls()
        obj.raw = line
        line = line.rstrip("\r\n")
        obj["key"], obj["value"] = line[2:].split("=", 1)
        return obj

    def __str__(self):
        return f"##{self['key']}={self['value']}"

    @staticmethod
    def is_type(raw: str) -> bool:
        if not raw.startswith("##"):
            return False
        key = raw[2:].split("=", 1)[0]
        return key not in ("INFO", "FILTER", "FORMAT", "contig")


class Fields(list):
    """The fields/column names"""

    kind = "fields"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.raw = None

    @classmethod
    def from_str(cls, line: str):
        obj = cls()
        obj.raw = line
        line = line.rstrip("\r\n")
        obj.extend(line[1:].split("\t"))
        return obj

    def __str__(self):
        return "#" + "\t".join(self)

    @property
    def samples(self):
        return self[9:]

    @staticmethod
    def is_type(raw: str) -> bool:
        return raw.startswith("#CHROM")


class Info(dict):
    """The INFO of the variant"""
    @classmethod
    def from_str(cls, infostr: str):
        obj = cls()
        for part in infostr.split(";"):
            # a flag
            if "=" not in part:
                obj[part] = True
            else:
                name, value = part.split("=", 1)
                obj[name] = value

        return obj

    def __str__(self) -> str:
        return ";".join(
            k if v is True else f"{k}={v}"
            for k, v in self.items()
            if v is not False
        )


class Format(list):
    """The FORMAT of the variant"""

    @classmethod
    def from_str(cls, formatstr: str):
        return cls(formatstr.split(":"))

    def __str__(self) -> str:
        return ":".join(self)


class Alt(list):
    """The ALT of the variant"""

    @classmethod
    def from_str(cls, altstr):
        return cls(altstr.split(","))

    def __str__(self) -> str:
        return ",".join(self)


class Filter(list):
    """The FILTER of the variant"""

    @classmethod
    def from_str(cls, filtstr: str):
        return cls(filtstr.split(";"))

    def __str__(self) -> str:
        return ";".join(self)


class Sample(dict):
    """One sample of the variant"""
    def __init__(self, values: Sequence[str], format: Format):
        super().__init__()
        self._format = format
        for name, value in zip(format, values):
            self[name] = value

    @property
    def format(self):
        return self._format

    @classmethod
    def from_str(cls, value_str: str, format: Format):
        return cls(value_str.split(":"), format)

    @classmethod
    def from_strs(cls, value_strs: Sequence[str], format: Format):
        return cls(value_strs, format)

    def __str__(self) -> str:
        values = [self[fmt] for fmt in self._format]
        return ":".join(values)


class Samples(list):
    """The samples of the variant"""

    def __init__(self, samples: Sequence[Sample], format: Format):
        super().__init__(samples)
        self._format = format

    @property
    def format(self):
        return self._format

    @classmethod
    def from_str(cls, sample_str: str, format: Format):
        return cls(
            [
                Sample.from_str(sam_str, format)
                for sam_str in sample_str.split("\t")
            ],
            format,
        )

    @classmethod
    def from_strs(cls, sample_strs: Sequence[str], format: Format):
        return cls(
            [
                Sample.from_str(sam_str, format)
                for sam_str in sample_strs
            ],
            format,
        )

    @classmethod
    def from_strss(cls, sample_strss: Sequence[Sequence[str]], format: Format):
        return cls(
            [
                Sample.from_strs(sam_strs, format)
                for sam_strs in sample_strss
            ],
            format,
        )

    def __str__(self) -> str:
        return "\t".join(str(s) for s in self)


class Variant:

    kind = "variant"

    def __init__(
        self,
        chrom: str,
        pos: int,
        id: str,
        ref: str,
        alt: Alt,
        qual: str,
        filter: Filter,
        info: Info,
        format: Format,
        samples: Samples,
    ):
        self.chrom = chrom
        self.pos = pos
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format = format
        self.samples = samples
        self.raw = None

    @classmethod
    def from_strs(
        cls,
        chrom: str,
        pos: int | str,
        id: str,
        ref: str,
        alt: str | Sequence[str],
        qual: str,
        filter: str | Sequence[str],
        info: str | dict,
        format: str | Sequence[str],
        samples: str | Sequence[str] | Sequence[Sequence[str]],
    ):
        format = (
            Format.from_str(format)
            if isinstance(format, str)
            else Format(format)
        )

        if isinstance(samples, str):
            samples = Samples.from_str(samples, format)
        elif isinstance(samples[0], str):
            samples = Samples.from_strs(samples, format)  # type: ignore
        else:
            samples = Samples.from_strss(samples, format)

        obj = cls(
            chrom,
            int(pos),
            id,
            ref,
            Alt.from_str(alt) if isinstance(alt, str) else Alt(alt),
            qual,
            Filter.from_str(filter)
            if isinstance(filter, str)
            else Filter(filter),
            Info.from_str(info) if isinstance(info, str) else Info(info),
            format,
            samples,
        )
        return obj

    @classmethod
    def from_str(cls, variant_line: str):
        raw = variant_line
        variant_line = variant_line.rstrip("\r\n")
        items = variant_line.split("\t")
        chrom = items[0]
        pos = int(items[1])
        id = items[2]
        ref = items[3]
        alt = Alt.from_str(items[4])
        qual = items[5]
        filter = Filter.from_str(items[6])
        info = Info.from_str(items[7])
        format = Format.from_str(items[8])
        samples = Samples.from_strs(items[9:], format)
        obj = cls(
            chrom,
            pos,
            id,
            ref,
            alt,
            qual,
            filter,
            info,
            format,
            samples,
        )
        obj.raw = raw
        return obj

    def __str__(self):
        return (
            f"{self.chrom}\t{self.pos}\t{self.id}\t{self.ref}\t"
            f"{self.alt}\t{self.qual}\t{self.filter}\t{self.info}\t"
            f"{self.format}\t{self.samples}"
        )

    def __repr__(self):
        return f"Variant({self.chrom}, {self.pos}, {self.id})"

    @staticmethod
    def is_type(raw: str) -> bool:
        return not raw.startswith("#")
