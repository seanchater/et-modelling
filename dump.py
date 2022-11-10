

            if isinstance(x, xr.Dataset):
                if "source" in x.encoding.keys():
                    # fp = x.encoding["source"]
                    # ds = xr.open_dataset(fp)
                    # ds = ds.close()
                    try:
                        # fp = xr.open_dataset(fp)
                        x = x.close()
                        os.remove(x)
                    except Exception as e:
                        print(x, " failed to remove:\n", e)
                    try:
                        # ds.close()
                        # fp = xr.open_dataset(fp)
                        x.close()
                        os.remove(x)
                    except Exception as e:
                        print(x, " failed to execute\n", e)
            elif isinstance(x, str):
                if os.path.isfile(x):
                    # fp = x
                    # ds = xr.open_dataset(fp)
                    # ds = ds.close()
                    try:
                        # ds.close()
                        # fp = xr.open_dataset(x)
                        x = x.close()
                        os.remove(x)
                    except Exception as e:
                        print(x, ' str vers. failed\n', e)