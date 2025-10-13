"""
System dependency checker for required bioinformatics tools.
"""

import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import logging


def check_command_availability(command: str) -> bool:
    """
    Check if a command is available in the system PATH.
    
    Args:
        command: Command name to check
        
    Returns:
        True if command is available, False otherwise
    """
    return shutil.which(command) is not None


def get_command_version(command: str, version_flag: str = '--version') -> Optional[str]:
    """
    Get version string for a command.
    
    Args:
        command: Command name
        version_flag: Flag to get version (default: --version)
        
    Returns:
        Version string or None if command not found
    """
    try:
        result = subprocess.run(
            [command, version_flag], 
            capture_output=True, 
            text=True, 
            timeout=10
        )
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            return result.stderr.strip()
    except (subprocess.TimeoutExpired, FileNotFoundError, subprocess.SubprocessError):
        return None


def check_python_packages(packages: List[str]) -> Dict[str, bool]:
    """
    Check if Python packages are installed.
    
    Args:
        packages: List of package names to check
        
    Returns:
        Dictionary mapping package names to availability status
    """
    results = {}
    for package in packages:
        try:
            __import__(package)
            results[package] = True
        except ImportError:
            results[package] = False
    
    return results


def check_r_packages(packages: List[str]) -> Dict[str, bool]:
    """
    Check if R packages are installed.
    
    Args:
        packages: List of R package names to check
        
    Returns:
        Dictionary mapping package names to availability status
    """
    results = {}
    
    if not check_command_availability('Rscript'):
        return {pkg: False for pkg in packages}
    
    for package in packages:
        try:
            r_command = f'if (!require("{package}", quietly=TRUE)) quit(status=1)'
            result = subprocess.run(
                ['Rscript', '-e', r_command],
                capture_output=True,
                timeout=30
            )
            results[package] = (result.returncode == 0)
        except subprocess.TimeoutExpired:
            results[package] = False
    
    return results


def verify_tools() -> Dict[str, Dict[str, any]]:
    """
    Verify availability of all required tools and their versions.
    
    Returns:
        Dictionary with tool status and versions
    """
    tools_status = {}
    
    # External bioinformatics tools
    external_tools = {
        'cutadapt': '--version',
        'vsearch': '--version', 
        'fastp': '--version',
        'Rscript': '--version'
    }
    
    for tool, version_flag in external_tools.items():
        available = check_command_availability(tool)
        version = get_command_version(tool, version_flag) if available else None
        
        tools_status[tool] = {
            'available': available,
            'version': version,
            'path': shutil.which(tool) if available else None
        }
    
    # Python packages
    python_packages = [
        'numpy', 'pandas', 'scipy', 'scikit-learn', 'biopython',
        'pysam', 'matplotlib', 'seaborn', 'tqdm', 'click',
        'datasketch', 'numba', 'joblib'
    ]
    
    python_status = check_python_packages(python_packages)
    tools_status['python_packages'] = python_status
    
    # R packages
    r_packages = ['dada2', 'Biostrings', 'ShortRead']
    r_status = check_r_packages(r_packages)
    tools_status['r_packages'] = r_status
    
    return tools_status


def check_dependencies(logger: Optional[logging.Logger] = None) -> Tuple[bool, Dict]:
    """
    Check all dependencies and log results.
    
    Args:
        logger: Optional logger for output
        
    Returns:
        Tuple of (all_available, status_dict)
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    logger.info("Checking system dependencies...")
    
    status = verify_tools()
    
    # Check critical dependencies
    critical_tools = ['cutadapt', 'vsearch', 'fastp', 'Rscript']
    critical_python = ['numpy', 'pandas', 'biopython', 'scipy']
    critical_r = ['dada2']
    
    missing_tools = []
    
    # Check external tools
    for tool in critical_tools:
        if tool in status and not status[tool]['available']:
            missing_tools.append(tool)
            logger.error(f"Critical tool not found: {tool}")
        elif tool in status and status[tool]['available']:
            logger.info(f"Found {tool}: {status[tool]['version']}")
    
    # Check Python packages
    python_status = status.get('python_packages', {})
    for pkg in critical_python:
        if not python_status.get(pkg, False):
            missing_tools.append(f"python:{pkg}")
            logger.error(f"Critical Python package not found: {pkg}")
        else:
            logger.info(f"Found Python package: {pkg}")
    
    # Check R packages
    r_status = status.get('r_packages', {})
    for pkg in critical_r:
        if not r_status.get(pkg, False):
            missing_tools.append(f"r:{pkg}")
            logger.error(f"Critical R package not found: {pkg}")
        else:
            logger.info(f"Found R package: {pkg}")
    
    all_available = len(missing_tools) == 0
    
    if all_available:
        logger.info("All critical dependencies are available!")
    else:
        logger.error(f"Missing dependencies: {missing_tools}")
    
    return all_available, status


def generate_dependency_report(output_path: Path) -> None:
    """
    Generate a detailed dependency report.
    
    Args:
        output_path: Path to save the report
    """
    status = verify_tools()
    
    with open(output_path, 'w') as f:
        f.write("# Dependency Check Report\\n\\n")
        f.write(f"Generated on: {Path(__file__).stat().st_mtime}\\n\\n")
        
        f.write("## External Tools\\n\\n")
        for tool, info in status.items():
            if isinstance(info, dict) and 'available' in info:
                status_icon = "✅" if info['available'] else "❌"
                f.write(f"- {status_icon} **{tool}**")
                if info['available']:
                    f.write(f": {info['version']} ({info['path']})")
                f.write("\\n")
        
        f.write("\\n## Python Packages\\n\\n")
        python_status = status.get('python_packages', {})
        for pkg, available in python_status.items():
            status_icon = "✅" if available else "❌"
            f.write(f"- {status_icon} {pkg}\\n")
        
        f.write("\\n## R Packages\\n\\n")
        r_status = status.get('r_packages', {})
        for pkg, available in r_status.items():
            status_icon = "✅" if available else "❌"
            f.write(f"- {status_icon} {pkg}\\n")