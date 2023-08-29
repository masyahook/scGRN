"""Download NDEx network from the web."""
import argparse
import json
import os

import ndex2

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download NDEx network using UUID.")
    parser.add_argument(
        "-i", "--id", type=str, help="The UUID of the NDEx network", required=True
    )
    parser.add_argument(
        "-f",
        "--folder",
        type=str,
        help="The folder where to save the network",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        help="The name of the .cx file, if None then the name from the NDEx "
        "will be used",
        default=None,
    )
    args = parser.parse_args()

    # Create NDEx2 python client
    client = ndex2.client.Ndex2()

    # Obtain metadata information about the network
    meta = client.get_network_summary(args.id)
    name = args.name if args.name is not None else meta["name"]

    # Obtain the network
    client_resp = client.get_network_as_cx_stream(args.id)

    # Convert downloaded network to NiceCXNetwork object
    net_cx = ndex2.create_nice_cx_from_raw_cx(json.loads(client_resp.content))

    # Save the network
    path_name = os.path.join(args.folder, f"{name}.cx")
    with open(path_name, "w") as f:
        json.dump(net_cx.to_cx(), f)

    print(f'Successfully saved the network at "{path_name}"!')
