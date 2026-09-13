"""Console entry points. Parse arguments before importing image analysis tools."""

import argparse


def alignment():
    """Run the existing alignment workflow from the active environment."""
    parser = argparse.ArgumentParser(description="Fibronectin alignment analysis")
    parser.add_argument("-i", "--input", required=True,
                        help="Path to a JSON file containing folder_paths")
    parser.add_argument("-a", "--angle_value", type=float, default=15,
                        help="Alignment angle in degrees (default: 15)")
    args = parser.parse_args()
    from .alignment_analysis import main_fibronectin_processing
    main_fibronectin_processing(args.input, args.angle_value)


def thickness():
    """Run the existing thickness workflow from the active environment."""
    parser = argparse.ArgumentParser(description="Fibronectin thickness analysis")
    parser.add_argument("-i", "--input", required=True,
                        help="Path to a JSON file containing folder_paths")
    args = parser.parse_args()
    from .thickness_analysis import main
    main(args.input)
