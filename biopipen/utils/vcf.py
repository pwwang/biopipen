from typing import Any, Sequence, Union


class HeaderItem(dict):

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

    kind = "contig"

    def __str__(self):
        return f"##contig=<ID={self['ID']}," f"length={self['length']}>"

    @staticmethod
    def is_type(raw: str) -> bool:
        return raw.startswith("##contig")


class HeaderGeneral(HeaderItem):

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

    @staticmethod
    def is_type(raw: str) -> bool:
        return raw.startswith("#CHROM")


class Info(dict):

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
            k
            if v is True
            else f"{k}={v}"
            for k, v in self.items()
            if v is not False
        )


class Format(list):

    @classmethod
    def from_str(cls, formatstr: str):
        return cls(formatstr.split(":"))

    def __str__(self) -> str:
        return ":".join(self)


class Alt(list):

    @classmethod
    def from_str(cls, altstr):
        return cls(altstr.split(","))

    def __str__(self) -> str:
        return ",".join(self)


class Filter(list):

    @classmethod
    def from_str(cls, filtstr: str):
        return cls(filtstr.split(";"))

    def __str__(self) -> str:
        return ";".join(self)


class Sample(Format):
    ...


class Samples(list):

    @classmethod
    def from_str(cls, sample_str: str):
        return cls(sample_str.split("\t"))

    @classmethod
    def from_strs(cls, sample_strs: list):
        return cls(Sample.from_str(s) for s in sample_strs)

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
        samples: Samples
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
        pos: Union[int, str],
        id: str,
        ref: str,
        alt: Union[str, Sequence[str]],
        qual: str,
        filter: Union[str, Sequence[str]],
        info: Union[str, dict],
        format: Union[str, Sequence[str]],
        samples: Union[str, Sequence[str], Sequence[Sequence[str]]],
    ):
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
            Info.from_str(info)
            if isinstance(info, str)
            else Info(info),
            Format.from_str(format)
            if isinstance(format, str)
            else Format(format),
            Samples.from_str(samples)
            if isinstance(samples, str)
            else Samples.from_strs(samples)
            if isinstance(samples, Sequence) and isinstance(samples[0], str)
            else Samples(samples),
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
        samples = Samples.from_strs(items[9:])
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
