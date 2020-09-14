"""Script entry point for bioprocs"""
from .args import params
from .modules import Module
from .commands import Command

def main():
    """Main entry point of biopipen CLI"""
    parsed = params.parse()
    klass = Command if Command.is_a(parsed.__command__) else Module
    klass.get(parsed.__command__).run(parsed[parsed.__command__])
