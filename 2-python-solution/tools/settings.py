# -*- coding: utf-8 -*-
"""Where this copy of the repository finds the public reference material, and how it runs the compiled program.

Everything that belongs to one machine is in one file at the repository's top, settings.cfg, which is not committed
(settings.cfg.example is; copy it and fill it in: 1-setup-and-documentation/README.md). Lines are "key = value"; a # starts a
comment. Keys:
  xpsi         the folder of X-PSI's source tree at commit bb5f461e605b22c6251b155b27a914a558757f93; only its
               examples/examples_fast/Data folder is read
  apjlab5968   the folder of the supplementary data of Bogdanov et al. 2019 (ApJ Letters 887, L26), 24 text files
  wsl          Windows only: the Linux distribution under WSL in which the compiled solution is built and run; left out,
               WSL's default distribution is used. On Linux the compiled program is run directly and the key is not read.
A Windows path (C:\\a\\b) read on the Linux side of a Windows machine becomes /mnt/c/a/b, and the other way round, so one
settings file serves both sides.
"""
import os
import re
import shlex
import subprocess

TOP = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SETTINGS_FILE = os.path.join(TOP, "settings.cfg")
DATA_SUBFOLDER = os.path.join("examples", "examples_fast", "Data")
HELP = "see 1-setup-and-documentation/README.md, section 'Tell the programs where the reference material is'"


class SettingsError(Exception):
    """A missing or wrong settings.cfg, or missing reference material. The message says what to do."""


def read_settings():
    if not os.path.isfile(SETTINGS_FILE):
        raise SettingsError("there is no settings file %s: copy settings.cfg.example to settings.cfg at the repository's top and "
                            "fill in the two folders (%s)" % (SETTINGS_FILE, HELP))
    out = {}
    with open(SETTINGS_FILE, encoding="utf-8-sig") as fh:          # -sig: a byte-order mark from a Windows editor is dropped
        for number, line in enumerate(fh, 1):
            text = line.split("#", 1)[0].strip()
            if not text:
                continue
            if "=" not in text:
                raise SettingsError("settings.cfg, line %d: expected 'key = value', found %r" % (number, text))
            key, value = text.split("=", 1)
            out[key.strip().lower()] = value.strip().strip('"').strip("'")
    return out


def this_side_path(path):
    """The path as the side that runs this Python sees it (Windows or the Linux side of a Windows machine)."""
    if os.name == "posix":
        m = re.match(r"^([A-Za-z]):[\\/](.*)$", path)
        if m:
            return "/mnt/%s/%s" % (m.group(1).lower(), m.group(2).replace("\\", "/"))
    else:
        m = re.match(r"^/mnt/([A-Za-z])/(.*)$", path)
        if m:
            return "%s:/%s" % (m.group(1).upper(), m.group(2))
    return path


def linux_side_path(path):
    """The path as the Linux side sees it (used for the commands sent to WSL from Windows)."""
    path = os.path.abspath(path)
    m = re.match(r"^([A-Za-z]):[\\/](.*)$", path)
    if m:
        return "/mnt/%s/%s" % (m.group(1).lower(), m.group(2).replace("\\", "/"))
    return path


def folder(key, must_hold):
    cfg = read_settings()
    if key not in cfg or not cfg[key]:
        raise SettingsError("settings.cfg has no line '%s = <folder>' (%s)" % (key, HELP))
    path = this_side_path(os.path.expanduser(cfg[key]))
    if not os.path.isdir(path):
        raise SettingsError("settings.cfg: %s = %s is not a folder on this side of the machine; check the spelling, and on Windows "
                            "write it as C:\\... or /mnt/c/... (%s)" % (key, path, HELP))
    if not os.path.isfile(os.path.join(path, must_hold)):
        raise SettingsError("settings.cfg: the folder %s = %s does not hold %s (%s)" % (key, path, must_hold, HELP))
    return path


def data_dir():
    """The public quick example's Data folder: X-PSI's examples/examples_fast/Data."""
    return os.path.join(folder("xpsi", os.path.join(DATA_SUBFOLDER, "xpsi_good_realisation.dat")), DATA_SUBFOLDER)


def ref_dir():
    """The folder of the agreed reference curves and tables (apjlab5968)."""
    return folder("apjlab5968", "OS1e_test_IM.txt")


def compiled_dir():
    return os.path.join(TOP, "3-compiled-solution")


def run_compiled(program, args, timeout=900):
    """Run one of the compiled programs (nicer_fast or nicer_fast_strict) in 3-compiled-solution with the given arguments and
    return what it printed. On Linux it runs directly; on Windows it runs in the WSL distribution named in settings.cfg (or
    WSL's default one). An empty answer means the program is not built or did not run."""
    cdir = compiled_dir()
    words = [str(a) for a in args]
    if os.name == "posix":
        exe = os.path.join(cdir, program)
        if not os.path.isfile(exe):
            return ""
        done = subprocess.run([exe] + words, cwd=cdir, capture_output=True, text=True, timeout=timeout)
        return done.stdout
    distro = read_settings().get("wsl", "")
    command = "cd %s && ./%s %s" % (shlex.quote(linux_side_path(cdir)), program, " ".join(shlex.quote(w) for w in words))
    argv = ["wsl.exe"] + (["-d", distro] if distro else []) + ["--", "bash", "-lc", command]
    done = subprocess.run(argv, capture_output=True, text=True, timeout=timeout)
    return done.stdout.replace("\x00", "")        # wsl.exe writes its output two bytes per character when it is not a console
